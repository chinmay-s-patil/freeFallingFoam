#!/usr/bin/env python3
"""
animate_motion.py — 6DoF Propeller Tumbling Animation in PyVista (GPU-Accelerated & Multi-Core Parallel)

Renders continuous or step-by-step 3D animations of the tumbling propeller
using PyVista, OpenGL hardware acceleration, and multi-process parallel rendering.

Playback Modes:
  1. Real-Time (--real-time):
     Treats 1 second of simulation time as exactly 1 second in the video.
  2. Time-Scaled (--video-sec-per-sim-sec X):
     Define how many video seconds correspond to 1 second of simulation time
     (e.g. 2.0 = 2x slow-mo, 0.5 = 2x fast-forward).
  3. Step-by-Step Raw Frames (--every-input):
     Treats every recorded simulation step in motion.dat as 1 video frame.
     Video length is directly controlled via --fps (Duration = Total Steps / FPS).

Parallelization:
  • Multi-worker GPU chunk rendering (-j / --workers N, default: 4).
  • Splits the frame sequence across N worker processes rendering concurrently.
  • Lossless, zero-overhead FFmpeg stream concatenation (-c copy) in < 0.2s.

Usage Examples:
  # 1. Real-time video with 4 parallel workers (default):
  python3 animate_motion.py --real-time

  # 2. Maximum speed on 8 parallel workers:
  python3 animate_motion.py --real-time -j 8

  # 3. Custom time scale (1s sim = 3s video, 3x slow-mo):
  python3 animate_motion.py --video-sec-per-sim-sec 3.0 -j 4

  # 4. Step-by-step raw frames at 120 FPS across 4 workers:
  python3 animate_motion.py --every-input --fps 120 -j 4

  # 5. Quick step-by-step preview (first 500 frames):
  python3 animate_motion.py --every-input --fps 60 --max-frames 500

  # 6. Single still screenshot:
  python3 animate_motion.py --preview
"""

import os
import sys
import argparse
import time
import shutil
import tempfile
import subprocess
import multiprocessing as mp
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation as R, Slerp

# Suppress harmless Matplotlib / 3D warnings
import warnings
warnings.filterwarnings("ignore")

import pyvista as pv


def locate_propeller_stl(case_dir):
    """Search common case and meshing directories for the propeller STL geometry."""
    candidates = [
        os.path.join(case_dir, "../../meshing/rotatingMesh/constant/triSurface/patch_2_1.stl"),
        os.path.join(case_dir, "../../meshing/rotatingMesh/constant/triSurface/MeshRefinement.stl"),
        os.path.join(case_dir, "constant/triSurface/propeller.stl"),
        os.path.join(case_dir, "constant/triSurface/patch_2_1.stl"),
    ]
    for c in candidates:
        abs_p = os.path.abspath(c)
        if os.path.exists(abs_p):
            return abs_p
    return None


def locate_motion_data(case_dir):
    """Locate tumbleMotion/motion.dat file."""
    p = os.path.join(case_dir, "postProcessing/tumbleMotion/motion.dat")
    if os.path.exists(p):
        return os.path.abspath(p)
    alt = os.path.join(case_dir, "processor0/postProcessing/tumbleMotion/motion.dat")
    if os.path.exists(alt):
        return os.path.abspath(alt)
    return None


def load_and_filter_motion(motion_path, skip_points=5):
    """
    Load attitude quaternions and time history from tumbleMotion/motion.dat.
    Filters restart duplicate timestamps and startup transients.
    """
    if not os.path.exists(motion_path):
        raise FileNotFoundError(f"Motion file not found: {motion_path}")

    # Read header and data
    with open(motion_path, "r") as f:
        first_line = f.readline().strip()

    if first_line.startswith("#"):
        raw_headers = [h for h in first_line.lstrip("#").strip().split() if h]
        headers = [h.split("(")[0].strip() for h in raw_headers]
        df = pd.read_csv(motion_path, sep=r"\s+", comment="#", names=headers)
    else:
        df = pd.read_csv(motion_path, sep=r"\s+")

    # Standardize column names
    col_map = {c: c.split("(")[0].replace("#", "").strip() for c in df.columns}
    df.rename(columns=col_map, inplace=True)

    times = df["Time"].values
    n = len(times)

    # Detect rewinds / duplicates
    diffs = np.diff(times)
    rewinds = np.where(diffs <= 0)[0]

    drop_indices = set()
    # Filter initial start
    if skip_points > 0 and n > 0:
        for k in range(min(skip_points, n)):
            drop_indices.add(df.index[k])

    for r_idx in rewinds:
        if times[r_idx] == times[r_idx + 1]:
            drop_indices.add(df.index[r_idx])
        for k in range(r_idx + 1, min(r_idx + 1 + skip_points, n)):
            drop_indices.add(df.index[k])

    df_clean = df.drop(index=list(drop_indices)).reset_index(drop=True)
    return df_clean


def build_slerp_interpolator(df_motion):
    """Construct a continuous Slerp interpolator from motion quaternion telemetry."""
    times = df_motion["Time"].values

    # Ensure strictly monotonic increasing time for Slerp
    diffs = np.diff(times, prepend=times[0] - 1.0)
    strict_idx = np.where(diffs > 1e-9)[0]
    times_clean = times[strict_idx]

    qw = df_motion["Quat_w"].values[strict_idx]
    qx = df_motion["Quat_x"].values[strict_idx]
    qy = df_motion["Quat_y"].values[strict_idx]
    qz = df_motion["Quat_z"].values[strict_idx]

    # SciPy Rotation expects [qx, qy, qz, qw]
    quats_xyzw = np.column_stack([qx, qy, qz, qw])
    norms = np.linalg.norm(quats_xyzw, axis=1, keepdims=True)
    quats_xyzw /= np.maximum(norms, 1e-12)

    # Enforce continuous sign convention: q and -q represent same orientation
    for i in range(1, len(quats_xyzw)):
        if np.dot(quats_xyzw[i], quats_xyzw[i - 1]) < 0.0:
            quats_xyzw[i] = -quats_xyzw[i]

    rotations = R.from_quat(quats_xyzw)
    slerp = Slerp(times_clean, rotations)

    # Telemetry for HUD
    rpm_mag = np.sqrt(df_motion["RPM_x"]**2 + df_motion["RPM_y"]**2 + df_motion["RPM_z"]**2).values[strict_idx]
    qw_clean = np.clip(np.abs(qw), -1.0, 1.0)
    tilt_3d = 2.0 * np.arccos(qw_clean) * (180.0 / np.pi)

    return slerp, times_clean, rpm_mag, tilt_3d


def extract_raw_step_telemetry(df_motion, stride=1, max_frames=None):
    """Extract raw step-by-step telemetry without interpolation (1 step = 1 frame)."""
    df_sub = df_motion.iloc[::stride]
    if max_frames is not None and max_frames > 0:
        df_sub = df_sub.iloc[:max_frames]

    times = df_sub["Time"].values
    qw = df_sub["Quat_w"].values
    qx = df_sub["Quat_x"].values
    qy = df_sub["Quat_y"].values
    qz = df_sub["Quat_z"].values

    quats_xyzw = np.column_stack([qx, qy, qz, qw])
    norms = np.linalg.norm(quats_xyzw, axis=1, keepdims=True)
    quats_xyzw /= np.maximum(norms, 1e-12)

    for i in range(1, len(quats_xyzw)):
        if np.dot(quats_xyzw[i], quats_xyzw[i - 1]) < 0.0:
            quats_xyzw[i] = -quats_xyzw[i]

    rotations = R.from_quat(quats_xyzw)
    rot_matrices = rotations.as_matrix()

    rpm_mag = np.sqrt(df_sub["RPM_x"]**2 + df_sub["RPM_y"]**2 + df_sub["RPM_z"]**2).values
    qw_clean = np.clip(np.abs(qw), -1.0, 1.0)
    tilt_3d = 2.0 * np.arccos(qw_clean) * (180.0 / np.pi)

    return times, rot_matrices, rpm_mag, tilt_3d


def _render_chunk_worker(worker_id, start_idx, end_idx, total_frames,
                         stl_path, eval_rot_matrices, eval_times, eval_tilt, eval_rpm,
                         fps, style, every_input, real_time, video_sec_per_sim_sec,
                         width, height, out_chunk_path, shared_counter, counter_lock):
    """
    Subprocess worker: renders an allocated slice of frames to a video part file.
    """
    try:
        # Load local mesh and feature edges
        mesh = pv.read(stl_path)
        edges = mesh.extract_feature_edges(feature_angle=25, boundary_edges=True, feature_edges=True, manifold_edges=False)
        base_points_mesh = mesh.points.copy()
        base_points_edges = edges.points.copy()

        pv.OFF_SCREEN = True
        plotter = pv.Plotter(off_screen=True, window_size=[width, height])
        plotter.set_background("#0d1117")

        plotter.camera_position = [(0.28, 0.20, 0.28), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
        plotter.camera.zoom(1.15)
        plotter.add_axes(interactive=False, line_width=2.0)

        mesh_actor = None
        if style in ["solid", "both"]:
            mesh_actor = plotter.add_mesh(
                mesh.copy(),
                color="#388bfd",
                opacity=1.0,  # opaque for maximum OpenGL throughput
                smooth_shading=True,
                specular=0.45,
                ambient=0.30
            )

        edges_actor = plotter.add_mesh(
            edges.copy(),
            color="#ffffff",
            line_width=2.5,
            render_lines_as_tubes=False
        )

        hud_actor = plotter.add_text(
            "INITIALIZING TELEMETRY...",
            name="hud_text",
            position="upper_left",
            font_size=11,
            color="#c9d1d9",
            font="courier"
        )

        plotter.open_movie(out_chunk_path, framerate=fps, quality=8)

        for idx in range(start_idx, end_idx):
            rot_mat = eval_rot_matrices[idx]
            t_sim = eval_times[idx]
            tilt_val = eval_tilt[idx]
            rpm_val = eval_rpm[idx]

            if mesh_actor is not None:
                mesh_actor.mapper.dataset.points = base_points_mesh @ rot_mat.T
            edges_actor.mapper.dataset.points = base_points_edges @ rot_mat.T

            if every_input:
                hud_text = f"STEP: {idx+1:,}/{total_frames:,} | SIM TIME: {t_sim:.5f}s | TILT: {tilt_val:5.1f}° | SPEED: {rpm_val:6.1f} RPM ({fps} FPS)"
            else:
                tag = "REAL-TIME" if real_time else f"{video_sec_per_sim_sec:.2f}x"
                hud_text = f"SIM TIME: {t_sim:.4f}s ({tag}) | TILT: {tilt_val:5.1f}° | SPEED: {rpm_val:6.1f} RPM ({fps} FPS)"

            hud_actor.set_text(2, hud_text)
            plotter.write_frame()

            if shared_counter is not None and counter_lock is not None:
                with counter_lock:
                    shared_counter.value += 1

        plotter.close()

    except Exception as e:
        print(f"\n[Worker {worker_id} Error] {e}", file=sys.stderr)
        raise e


def render_animation(stl_path, df_motion, output_path="propeller_tumble_animation.mp4",
                     fps=60, real_time=False, video_sec_per_sim_sec=1.0,
                     every_input=False, stride=1, max_frames=None, workers=4,
                     style="solid", preview=False, width=1280, height=720):
    """
    Render 3D animated propeller video or still preview using PyVista on the GPU,
    parallelized across multiple worker processes.
    """
    print("=" * 72)
    print("  PyVista GPU 6DoF Propeller Tumbling Animation (Multi-Core Parallel)")
    print("=" * 72)

    # 1. GPU Information
    gpu_info = pv.GPUInfo()
    print(f"  • Hardware GPU    : {gpu_info.renderer} ({gpu_info.version})")
    print(f"  • Geometry STL    : {os.path.basename(stl_path)}")
    print(f"  • Telemetry Rows  : {len(df_motion):,} total filtered steps")

    # 2. Determine Frame Sequence & Timing
    if every_input:
        eval_times, eval_rot_matrices, eval_rpm, eval_tilt = extract_raw_step_telemetry(
            df_motion, stride=stride, max_frames=max_frames
        )
        total_frames = len(eval_times)
        duration_sim = eval_times[-1] - eval_times[0] if total_frames > 1 else 0.0
        duration_video = total_frames / float(fps)
        mode_label = f"STEP-BY-STEP (--every-input, 1 input = 1 frame, stride={stride})"

        print(f"  • Playback Mode   : {mode_label}")
        print(f"  • Simulation Time : {eval_times[0]:.4f}s to {eval_times[-1]:.4f}s (ΔT = {duration_sim:.4f}s)")
        print(f"  • Target FPS      : {fps} FPS")
        print(f"  • Total Frames    : {total_frames:,} frames")
        print(f"  • Video Duration  : {duration_video:.2f} s ({duration_video/60:.2f} minutes)")

    else:
        if real_time:
            time_scale = 1.0
            mode_label = "REAL-TIME (1s sim time = 1s video)"
        else:
            time_scale = float(video_sec_per_sim_sec)
            mode_label = f"TIME-SCALED (1s sim time = {time_scale:.2f}s video, {time_scale:.2f}x speed factor)"

        slerp, times_clean, rpm_series, tilt_series = build_slerp_interpolator(df_motion)
        t_start, t_end = times_clean[0], times_clean[-1]
        duration_sim = t_end - t_start
        duration_video = max(duration_sim * time_scale, 1.0 / fps)
        total_frames = max(1, int(np.round(duration_video * fps)))

        if max_frames is not None and max_frames > 0:
            total_frames = min(total_frames, max_frames)
            duration_video = total_frames / float(fps)
            t_span_eval = total_frames / (fps * time_scale)
            eval_times = np.linspace(t_start, min(t_start + t_span_eval, t_end), total_frames)
        else:
            eval_times = np.linspace(t_start, t_end, total_frames)

        eval_rotations = slerp(eval_times)
        eval_rot_matrices = eval_rotations.as_matrix()
        eval_rpm = np.interp(eval_times, times_clean, rpm_series)
        eval_tilt = np.interp(eval_times, times_clean, tilt_series)

        print(f"  • Playback Mode   : {mode_label}")
        print(f"  • Simulation Time : {t_start:.4f}s to {t_end:.4f}s (ΔT = {duration_sim:.4f}s)")
        print(f"  • Target FPS      : {fps} FPS")
        print(f"  • Total Frames    : {total_frames:,} frames")
        print(f"  • Video Duration  : {duration_video:.2f} seconds")

    print(f"  • Render Style    : {style.upper()} (Cobalt blue with glowing white outline)")

    # ----------------------------------------------------
    # PREVIEW MODE: Save single high-res still image
    # ----------------------------------------------------
    if preview:
        mesh = pv.read(stl_path)
        edges = mesh.extract_feature_edges(feature_angle=25, boundary_edges=True, feature_edges=True, manifold_edges=False)
        base_points_mesh = mesh.points.copy()
        base_points_edges = edges.points.copy()

        pv.OFF_SCREEN = True
        plotter = pv.Plotter(off_screen=True, window_size=[width, height])
        plotter.set_background("#0d1117")
        plotter.camera_position = [(0.28, 0.20, 0.28), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
        plotter.camera.zoom(1.15)
        plotter.add_axes(interactive=False, line_width=2.0)

        mesh_actor = None
        if style in ["solid", "both"]:
            mesh_actor = plotter.add_mesh(mesh.copy(), color="#388bfd", opacity=1.0, smooth_shading=True, specular=0.45, ambient=0.30)
        edges_actor = plotter.add_mesh(edges.copy(), color="#ffffff", line_width=2.5)

        hud_actor = plotter.add_text("INIT", position="upper_left", font_size=11, color="#c9d1d9", font="courier")

        mid_idx = total_frames // 2
        rot_mat = eval_rot_matrices[mid_idx]
        if mesh_actor is not None:
            mesh_actor.mapper.dataset.points = base_points_mesh @ rot_mat.T
        edges_actor.mapper.dataset.points = base_points_edges @ rot_mat.T

        if every_input:
            hud_text = f"STEP: {mid_idx+1:,}/{total_frames:,} | SIM TIME: {eval_times[mid_idx]:.5f}s | TILT: {eval_tilt[mid_idx]:.1f}° | SPEED: {eval_rpm[mid_idx]:.1f} RPM ({fps} FPS)"
        else:
            tag = "REAL-TIME" if real_time else f"{video_sec_per_sim_sec:.2f}x"
            hud_text = f"SIM TIME: {eval_times[mid_idx]:.4f}s ({tag}) | TILT: {eval_tilt[mid_idx]:.1f}° | SPEED: {eval_rpm[mid_idx]:.1f} RPM ({fps} FPS)"

        hud_actor.set_text(2, hud_text)
        preview_path = "propeller_tumble_preview.png"
        plotter.screenshot(preview_path)
        plotter.close()
        print(f"\n[Success] Still preview image rendered to: {os.path.abspath(preview_path)}")
        return

    # Determine worker count
    num_workers = max(1, int(workers))
    # Avoid spawning more workers than frames
    num_workers = min(num_workers, total_frames)

    # ----------------------------------------------------
    # SEQUENTIAL EXECUTION (workers <= 1)
    # ----------------------------------------------------
    if num_workers <= 1:
        print(f"  • Workers         : 1 (Sequential mode)")
        print(f"\nEncoding video to: {output_path} ...")
        t0_render = time.time()
        _render_chunk_worker(
            worker_id=0,
            start_idx=0,
            end_idx=total_frames,
            total_frames=total_frames,
            stl_path=stl_path,
            eval_rot_matrices=eval_rot_matrices,
            eval_times=eval_times,
            eval_tilt=eval_tilt,
            eval_rpm=eval_rpm,
            fps=fps,
            style=style,
            every_input=every_input,
            real_time=real_time,
            video_sec_per_sim_sec=video_sec_per_sim_sec,
            width=width,
            height=height,
            out_chunk_path=output_path,
            shared_counter=None,
            counter_lock=None
        )
        dt_total = time.time() - t0_render
        print(f"\n[Success] Animation rendered in {dt_total:.1f}s (Average: {total_frames/max(dt_total,1e-3):.1f} FPS)!")
        print(f"Video file saved: {os.path.abspath(output_path)} ({os.path.getsize(output_path)/1024/1024:.2f} MB)")
        return

    # ----------------------------------------------------
    # MULTI-PROCESS PARALLEL EXECUTION (workers > 1)
    # ----------------------------------------------------
    print(f"  • Parallel Workers: {num_workers} processes (Concurrent GPU rendering)")
    temp_dir = tempfile.mkdtemp(prefix="pv_tumble_parts_")

    try:
        # Split frame indices into balanced contiguous chunks
        splits = np.array_split(np.arange(total_frames), num_workers)
        chunk_tasks = []
        chunk_files = []

        for cid, c_range in enumerate(splits):
            if len(c_range) == 0:
                continue
            chunk_file = os.path.join(temp_dir, f"part_{cid:03d}.mp4")
            chunk_files.append(chunk_file)
            chunk_tasks.append((cid, c_range[0], c_range[-1] + 1, chunk_file))

        active_workers = len(chunk_tasks)
        print(f"\nDispatching {total_frames:,} frames across {active_workers} parallel workers ...")

        shared_counter = mp.Value('i', 0)
        counter_lock = mp.Lock()

        processes = []
        for (cid, s_idx, e_idx, c_file) in chunk_tasks:
            p = mp.Process(
                target=_render_chunk_worker,
                args=(
                    cid, s_idx, e_idx, total_frames,
                    stl_path, eval_rot_matrices, eval_times, eval_tilt, eval_rpm,
                    fps, style, every_input, real_time, video_sec_per_sim_sec,
                    width, height, c_file, shared_counter, counter_lock
                )
            )
            processes.append(p)

        t0_render = time.time()
        for p in processes:
            p.start()

        # Monitor progress across all workers
        last_print = 0
        while any(p.is_alive() for p in processes):
            now = time.time()
            with counter_lock:
                current = shared_counter.value

            if now - last_print >= 0.5 or current == total_frames:
                pct = (current / total_frames) * 100.0
                elapsed = now - t0_render
                fps_render = current / max(elapsed, 1e-3)
                eta = (total_frames - current) / max(fps_render, 1e-3)
                sys.stdout.write(
                    f"\r  Rendered {current:,}/{total_frames:,} frames ({pct:5.1f}%) | "
                    f"Speed: {fps_render:5.1f} FPS | ETA: {eta:4.1f}s [{active_workers} workers]"
                )
                sys.stdout.flush()
                last_print = now
            time.sleep(0.1)

        for p in processes:
            p.join()
            if p.exitcode != 0:
                raise RuntimeError(f"Worker process failed with exit code {p.exitcode}")

        dt_render = time.time() - t0_render
        print(f"\n\nAll {active_workers} chunks completed in {dt_render:.1f}s ({total_frames/max(dt_render,1e-3):.1f} FPS overall)!")

        # ----------------------------------------------------
        # Stitch chunks with FFmpeg stream copy (Instant & Lossless)
        # ----------------------------------------------------
        print("Stitching video chunks losslessly with FFmpeg stream copy ...")
        concat_list_file = os.path.join(temp_dir, "concat_list.txt")
        with open(concat_list_file, "w") as f:
            for c_path in chunk_files:
                f.write(f"file '{c_path}'\n")

        ffmpeg_cmd = [
            "ffmpeg", "-y", "-v", "error",
            "-f", "concat", "-safe", "0",
            "-i", concat_list_file,
            "-c", "copy",
            output_path
        ]
        subprocess.run(ffmpeg_cmd, check=True)

        total_time = time.time() - t0_render
        print(f"[Success] Video saved: {os.path.abspath(output_path)} ({os.path.getsize(output_path)/1024/1024:.2f} MB)")
        print(f"Total time elapsed: {total_time:.1f}s (Effective speed: {total_frames/max(total_time,1e-3):.1f} FPS)")

    finally:
        # Clean up temporary part files
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(
        description="6DoF Propeller Tumbling Animation in PyVista (GPU Accelerated & Multi-Core Parallel)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 1. Real-time playback on 4 parallel workers (default):
  python3 animate_motion.py --real-time

  # 2. Maximum parallel speed on 8 CPU workers:
  python3 animate_motion.py --real-time -j 8

  # 3. Custom time scale (1s sim time = 3.5s video):
  python3 animate_motion.py --video-sec-per-sim-sec 3.5 -j 4

  # 4. Every simulation input step as 1 frame at 120 FPS:
  python3 animate_motion.py --every-input --fps 120 -j 6

  # 5. Step-by-step test preview (first 500 frames):
  python3 animate_motion.py --every-input --fps 60 --max-frames 500
"""
    )
    parser.add_argument("--case", type=str, default=".", help="Case directory (default: current directory)")
    parser.add_argument("-o", "--output", type=str, default="propeller_tumble_animation.mp4", help="Output MP4 filename")

    # Multi-Processing Worker Argument
    parser.add_argument(
        "-j", "--workers", "--parallel",
        type=int,
        default=4,
        dest="workers",
        help="Number of parallel worker processes for rendering (default: 4, set to 1 for sequential)"
    )

    # Timing Mode Arguments
    parser.add_argument(
        "--real-time",
        action="store_true",
        help="Treat 1 second of simulation time as 1 second in the video (real-time playback)"
    )
    parser.add_argument(
        "--video-sec-per-sim-sec", "--time-scale", "--sec-per-sec", "--slowdown",
        type=float,
        default=None,
        dest="video_sec_per_sim_sec",
        help="Set how many video seconds correspond to 1 second of simulation time (e.g. 1.0 = real-time, 2.0 = 2s video per 1s sim [2x slow-mo], 0.5 = 2x fast-forward). Default: 1.0"
    )
    parser.add_argument(
        "--every-input", "--raw-frames", "--every-step",
        action="store_true",
        dest="every_input",
        help="Treat every recorded simulation input step in motion.dat as 1 video frame (no interpolation). Total video duration is determined by total steps / FPS."
    )

    # Frame Rate and Sub-sampling Arguments
    parser.add_argument(
        "--fps",
        type=int,
        default=60,
        help="Output video framerate (default: 60 FPS). In --every-input mode, this directly controls playback speed and total video length."
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Step stride in --every-input mode (e.g. 1 = every step, 2 = every 2nd step, 10 = every 10th step). Default: 1"
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Maximum number of frames to render (useful for quick previews or testing)"
    )

    # Aesthetics & Filtering
    parser.add_argument(
        "--style",
        type=str,
        default="solid",
        choices=["solid", "outline"],
        help="Render style: 'solid' (colored body + sharp outline) or 'outline' (wireframe perimeter)"
    )
    parser.add_argument("--preview", action="store_true", help="Generate a still preview image instead of full video")
    parser.add_argument("--skip-restart-points", type=int, default=5, help="Number of restart transient points to filter (default: 5)")
    args = parser.parse_args()

    # Resolve timing parameters
    if args.real_time:
        video_sec_per_sim_sec = 1.0
    elif args.video_sec_per_sim_sec is not None:
        video_sec_per_sim_sec = args.video_sec_per_sim_sec
    else:
        video_sec_per_sim_sec = 1.0  # Default to real-time if not specified

    case_dir = os.path.abspath(args.case)
    stl_path = locate_propeller_stl(case_dir)
    motion_path = locate_motion_data(case_dir)

    if not stl_path:
        print("[Error] Could not locate propeller STL geometry.", file=sys.stderr)
        sys.exit(1)

    if not motion_path:
        print(f"[Error] Could not locate motion.dat in: {case_dir}", file=sys.stderr)
        sys.exit(1)

    df_motion = load_and_filter_motion(motion_path, skip_points=args.skip_restart_points)
    out_file = os.path.join(case_dir, args.output)

    render_animation(
        stl_path=stl_path,
        df_motion=df_motion,
        output_path=out_file,
        fps=args.fps,
        real_time=args.real_time,
        video_sec_per_sim_sec=video_sec_per_sim_sec,
        every_input=args.every_input,
        stride=args.stride,
        max_frames=args.max_frames,
        workers=args.workers,
        style=args.style,
        preview=args.preview
    )


if __name__ == "__main__":
    main()
