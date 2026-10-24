// ============================================================
//  predictRotorDynamics.H — 6DoF Rotational Predictor Step
// ============================================================

qBody_old = qBody;
omegaB_old = omegaB;

// Advance orientation by dt using omegaB from previous time step
// so that mesh.controlledUpdate() moves the inner mesh for PIMPLE iterations
const scalar dtPredict = runTime.deltaTValue();
const scalar wMagPredict = mag(omegaB);

if (wMagPredict > 1e-12 && dtPredict > 1e-15)
{
    const quaternion dq(omegaB/wMagPredict, wMagPredict*dtPredict);
    qBody = qBody_old * dq;
    qBody.normalise();
}

// Update uniform registered fields so tFFRotationalMotion reads the predicted orientation
attitudeW.value() = qBody.w();
attitudeVec.value() = qBody.v();
