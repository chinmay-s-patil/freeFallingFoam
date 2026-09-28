#!/usr/bin/env python3
"""
plot_motion.py
--------------
Python tracking script for freeFallingFoam (steelBallRun).
Parses postProcessing/fallingMotion/motion.dat and generates live updating graphs of:
1. Velocity vs Time
2. Acceleration vs Time
3. All Forces (Gravity, Aero, Pressure, Viscous, Net Force)

Features:
- Live popup GUI with automatic refresh (default every 30s)
- Interactive 'Refresh Now' button for manual updates
- Image export (motion_tracking.png)

Usage:
	python3 plot_motion.py [--file PATH] [--save FILENAME] [--interval SEC] [--no-popup]
"""

import sys
import os
import argparse
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

def parse_motion_dat(filepath):
	"""Parse postProcessing/fallingMotion/motion.dat into a structured numpy dict."""
	if not os.path.exists(filepath):
		raise FileNotFoundError(f"Data file not found at: {filepath}")

	with open(filepath, 'r') as f:
		lines = [line.strip() for line in f if line.strip()]

	if not lines:
		raise ValueError("motion.dat file is empty.")

	data_lines = []
	for line in lines:
		if line.startswith("#"):
			continue
		parts = line.split()
		# Filter for valid rows with expected 27 scalar fields
		if len(parts) == 27:
			try:
				row = [float(val) for val in parts]
				data_lines.append(row)
			except ValueError:
				continue

	if not data_lines:
		raise ValueError("No valid numeric data lines found in motion.dat.")

	data = np.array(data_lines)

	# Standard column mapping based on fallingMotion.C / updateFallingFrameMotion.H
	# Header: Time Altitude Vx Vy Vz Ax Ay Az g_x g_y g_z rho Fpress_x Fpress_y Fpress_z Fvisc_x Fvisc_y Fvisc_z Faero_x Faero_y Faero_z Fgrav_x Fgrav_y Fgrav_z Fnet_x Fnet_y Fnet_z
	cols = {
		'time': data[:, 0],
		'alt': data[:, 1],
		'Vx': data[:, 2], 'Vy': data[:, 3], 'Vz': data[:, 4],
		'Ax': data[:, 5], 'Ay': data[:, 6], 'Az': data[:, 7],
		'gx': data[:, 8], 'gy': data[:, 9], 'gz': data[:, 10],
		'rho': data[:, 11],
		'Fpress_x': data[:, 12], 'Fpress_y': data[:, 13], 'Fpress_z': data[:, 14],
		'Fvisc_x': data[:, 15], 'Fvisc_y': data[:, 16], 'Fvisc_z': data[:, 17],
		'Faero_x': data[:, 18], 'Faero_y': data[:, 19], 'Faero_z': data[:, 20],
		'Fgrav_x': data[:, 21], 'Fgrav_y': data[:, 22], 'Fgrav_z': data[:, 23],
		'Fnet_x': data[:, 24], 'Fnet_y': data[:, 25], 'Fnet_z': data[:, 26],
	}

	# Computed scalar magnitudes
	cols['V_mag'] = np.sqrt(cols['Vx']**2 + cols['Vy']**2 + cols['Vz']**2)
	cols['A_mag'] = np.sqrt(cols['Ax']**2 + cols['Ay']**2 + cols['Az']**2)
	cols['Fgrav_mag'] = np.sqrt(cols['Fgrav_x']**2 + cols['Fgrav_y']**2 + cols['Fgrav_z']**2)
	cols['Faero_mag'] = np.sqrt(cols['Faero_x']**2 + cols['Faero_y']**2 + cols['Faero_z']**2)
	cols['Fnet_mag'] = np.sqrt(cols['Fnet_x']**2 + cols['Fnet_y']**2 + cols['Fnet_z']**2)

	return cols


class LiveMotionViewer:
	def __init__(self, filepath, save_filename="motion_tracking.png", interval=30, show_popup=True):
		self.filepath = filepath
		self.save_filename = save_filename
		self.interval_sec = interval
		self.show_popup = show_popup

		plt.style.use('seaborn-v0_8-darkgrid' if 'seaborn-v0_8-darkgrid' in plt.style.available else 'default')

		self.fig, self.axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
		self.fig.subplots_adjust(bottom=0.11, top=0.93, hspace=0.35)

		# Add Button for Manual Refresh
		self.ax_btn = self.fig.add_axes([0.80, 0.02, 0.16, 0.045])
		self.btn_refresh = Button(self.ax_btn, ' Refresh Now ', color='#e1f5fe', hovercolor='#b3e5fc')
		self.btn_refresh.on_clicked(self.on_refresh_click)

		# Status text at bottom left
		self.status_text = self.fig.text(0.08, 0.03, "Initializing...", fontsize=10, color='#555555', verticalalignment='center')

		# Perform initial data load and plotting
		self.refresh_data()

		# Set up timer for auto-update if popup mode and interval > 0
		if self.show_popup and self.interval_sec > 0:
			self.timer = self.fig.canvas.new_timer(interval=self.interval_sec * 1000)
			self.timer.add_callback(self.refresh_data)
			self.timer.start()

	def on_refresh_click(self, event):
		self.refresh_data()

	def refresh_data(self):
		now_str = datetime.datetime.now().strftime("%H:%M:%S")
		try:
			cols = parse_motion_dat(self.filepath)
			self.draw_plots(cols)
			if self.save_filename:
				plt.savefig(self.save_filename, dpi=200)
				print(f"[{now_str}] Plot updated and saved to: {self.save_filename}")
			self.status_text.set_text(f"Last updated: {now_str} | Auto-refresh: {self.interval_sec}s")
			self.status_text.set_color('#2e7d32')
		except Exception as e:
			err_msg = f"Update failed ({now_str}): {e}"
			self.status_text.set_text(err_msg)
			self.status_text.set_color('#c62828')
			print(f"[WARNING] [{now_str}] {e}", file=sys.stderr)

		if self.show_popup:
			self.fig.canvas.draw_idle()

	def draw_plots(self, cols):
		ax_v, ax_a, ax_f = self.axes

		ax_v.clear()
		ax_a.clear()
		ax_f.clear()

		self.fig.suptitle("steelBallRun — Free Fall Dynamics Tracking", fontsize=15, fontweight='bold', y=0.97)

		t = cols['time']
		# Step-wise velocity difference: ΔV(t) = V(t) - V(t-1)
		delta_Vy = np.diff(cols['Vy'], prepend=cols['Vy'][0])
		curr_Vy = cols['Vy'][-1]
		curr_Vmag = cols['V_mag'][-1]

		# 1. Step Velocity Change (Delta V) Plot
		ax_v.plot(t, delta_Vy, color='#1f77b4', linewidth=1.8, linestyle='-', label='ΔVy = Vy(t) - Vy(t-1) [m/s]')
		ax_v.set_ylabel("ΔVy per step [m/s]", fontsize=11, fontweight='bold')
		ax_v.set_title("Step Velocity Change ΔV(t) = V(t) - V(t-1)", fontsize=12, fontweight='bold', loc='left')

		# Prominently show current falling velocity
		ax_v.text(0.02, 0.85, f"Current Falling Velocity:\nVy = {curr_Vy:.4f} m/s  (|V| = {curr_Vmag:.4f} m/s)",
		          transform=ax_v.transAxes, fontsize=10, fontweight='bold',
		          bbox=dict(boxstyle='round,pad=0.4', facecolor='#ffffff', edgecolor='#1f77b4', alpha=0.9),
		          verticalalignment='top')

		ax_v.legend(loc='upper right', frameon=True)
		ax_v.grid(True, linestyle=':', alpha=0.7)

		# 2. Acceleration Plot
		ax_a.plot(t, cols['Ay'], color='#ffbb78', linewidth=1.5, linestyle='--', label='Vertical Accel Ay (m/s²)')
		# ax_a.plot(t, cols['A_mag'], color='#ff7f0e', linewidth=1.8, linestyle='-', label='Accel Mag |a| (m/s²)')
		ax_a.set_ylabel("Acceleration [m/s²]", fontsize=11, fontweight='bold')
		ax_a.set_title("Acceleration Tracking", fontsize=12, fontweight='bold', loc='left')
		ax_a.legend(loc='upper right', frameon=True)
		ax_a.grid(True, linestyle=':', alpha=0.7)

		# 3. Forces Plot
		ax_f.plot(t, cols['Fgrav_mag'], color='#2ca02c', linewidth=1.8, linestyle='--', label='Gravity |Fgrav| [N]', alpha=0.85)
		ax_f.plot(t, cols['Faero_y'], color='#9467bd', linewidth=1.8, linestyle='-.', label='Aero Drag Faero [N]', alpha=0.85)
		ax_f.plot(t, cols['Fpress_y'], color='#17becf', linewidth=1.2, linestyle=':', label='Pressure Fpress [N]', alpha=0.7)
		# ax_f.plot(t, cols['Fvisc_y'], color='#8c564b', linewidth=1.2, linestyle=':', label='Viscous Fvisc [N]', alpha=0.7)
		# ax_f.plot(t, cols['Fnet_y'], color='#d62728', linewidth=2.5, linestyle='-', label='Net Force Fnet [N]', zorder=10)

		ax_f.set_xlabel("Time [s]", fontsize=11, fontweight='bold')
		ax_f.set_ylabel("Force [N]", fontsize=11, fontweight='bold')
		ax_f.set_title("Forces & Net Force Tracking", fontsize=12, fontweight='bold', loc='left')
		ax_f.legend(loc='upper right', frameon=True, facecolor='#ffffff', edgecolor='#d62728', framealpha=0.95)
		ax_f.grid(True, linestyle=':', alpha=0.7)
		# ax_f.set_yscale('symlog', linthresh=1e-2)

	def show(self):
		if self.show_popup:
			plt.show()

def main():
	default_file = "postProcessing/fallingMotion/motion.dat"
	if not os.path.exists(default_file) and os.path.exists("processor0/postProcessing/fallingMotion/motion.dat"):
		default_file = "processor0/postProcessing/fallingMotion/motion.dat"

	parser = argparse.ArgumentParser(description="Track steelBallRun free fall velocity, acceleration, and forces.")
	parser.add_argument("--file", default=default_file, help="Path to motion.dat")
	parser.add_argument("--save", default="motion_tracking.png", help="Output plot filename")
	parser.add_argument("--interval", type=int, default=30, help="Auto-refresh interval in seconds (default: 30)")
	parser.add_argument("--no-popup", action="store_true", help="Disable GUI popup window and only save image once")
	args = parser.parse_args()

	try:
		viewer = LiveMotionViewer(
			filepath=args.file,
			save_filename=args.save,
			interval=args.interval,
			show_popup=not args.no_popup
		)
		viewer.show()
	except Exception as e:
		print(f"[ERROR] Failed to process tracking data: {e}", file=sys.stderr)
		sys.exit(1)

if __name__ == "__main__":
	main()


