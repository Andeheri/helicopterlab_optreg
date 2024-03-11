/*
 * helicopter_dag_3.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "helicopter_dag_3".
 *
 * Model version              : 11.13
 * Simulink Coder version : 9.4 (R2020b) 29-Jul-2020
 * C source code generated on : Mon Mar 11 15:37:15 2024
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helicopter_dag_3.h"
#include "helicopter_dag_3_private.h"
#include "helicopter_dag_3_dt.h"

/* Block signals (default storage) */
B_helicopter_dag_3_T helicopter_dag_3_B;

/* Continuous states */
X_helicopter_dag_3_T helicopter_dag_3_X;

/* Block states (default storage) */
DW_helicopter_dag_3_T helicopter_dag_3_DW;

/* Real-time model */
static RT_MODEL_helicopter_dag_3_T helicopter_dag_3_M_;
RT_MODEL_helicopter_dag_3_T *const helicopter_dag_3_M = &helicopter_dag_3_M_;

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = (ODE3_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  helicopter_dag_3_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  helicopter_dag_3_output();
  helicopter_dag_3_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  helicopter_dag_3_output();
  helicopter_dag_3_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter_dag_3_output(void)
{
  /* local block i/o variables */
  real_T rtb_Frontgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_DeadZoney;
  real_T rtb_DeadZonex;
  real_T rtb_Backgain;
  real_T rtb_Gain1_idx_2;
  real_T rtb_Gain1_idx_3;
  real_T *lastU;
  int8_T rtAction;
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
    /* set solver stop time */
    if (!(helicopter_dag_3_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter_dag_3_M->solverInfo,
                            ((helicopter_dag_3_M->Timing.clockTickH0 + 1) *
        helicopter_dag_3_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter_dag_3_M->solverInfo,
                            ((helicopter_dag_3_M->Timing.clockTick0 + 1) *
        helicopter_dag_3_M->Timing.stepSize0 +
        helicopter_dag_3_M->Timing.clockTickH0 *
        helicopter_dag_3_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter_dag_3_M)) {
    helicopter_dag_3_M->Timing.t[0] = rtsiGetT(&helicopter_dag_3_M->solverInfo);
  }

  /* Reset subsysRan breadcrumbs */
  srClearBC(helicopter_dag_3_DW.IfActionSubsystem_SubsysRanBC);
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter_dag_3/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (helicopter_dag_3_DW.HILReadEncoderTimebase_Task, 1,
         &helicopter_dag_3_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter_dag_3_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_DeadZoney = helicopter_dag_3_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_DeadZonex = helicopter_dag_3_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helicopter_dag_3_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter_dag_3_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helicopter_dag_3_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helicopter_dag_3_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicopter_dag_3_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter_dag_3_B.FromWorkspace[0])[elIdx] =
                pDataValues[currTimeIndex];
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter_dag_3_B.FromWorkspace[0])[elIdx] =
                pDataValues[currTimeIndex + 1];
              pDataValues += 141;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 4; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&helicopter_dag_3_B.FromWorkspace[0])[elIdx] = (real_T)
              rtInterpolate(d1, d2, f1, f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
    /* Gain: '<S4>/Travel: Count to rad' incorporates:
     *  Gain: '<S4>/Travel_gain'
     */
    helicopter_dag_3_B.TravelCounttorad = helicopter_dag_3_P.travel_gain *
      rtb_HILReadEncoderTimebase_o1 * helicopter_dag_3_P.TravelCounttorad_Gain;

    /* Gain: '<S14>/Gain' */
    helicopter_dag_3_B.Gain = helicopter_dag_3_P.Gain_Gain *
      helicopter_dag_3_B.TravelCounttorad;

    /* Gain: '<S4>/Pitch: Count to rad' */
    helicopter_dag_3_B.PitchCounttorad = helicopter_dag_3_P.PitchCounttorad_Gain
      * rtb_DeadZoney;

    /* Gain: '<S11>/Gain' */
    helicopter_dag_3_B.Gain_i = helicopter_dag_3_P.Gain_Gain_a *
      helicopter_dag_3_B.PitchCounttorad;
  }

  /* Gain: '<S15>/Gain' incorporates:
   *  TransferFcn: '<S4>/Travel: Transfer Fcn'
   */
  helicopter_dag_3_B.Gain_d = (helicopter_dag_3_P.TravelTransferFcn_C *
    helicopter_dag_3_X.TravelTransferFcn_CSTATE +
    helicopter_dag_3_P.TravelTransferFcn_D * helicopter_dag_3_B.TravelCounttorad)
    * helicopter_dag_3_P.Gain_Gain_l;

  /* Gain: '<S12>/Gain' incorporates:
   *  TransferFcn: '<S4>/Pitch: Transfer Fcn'
   */
  helicopter_dag_3_B.Gain_b = (helicopter_dag_3_P.PitchTransferFcn_C *
    helicopter_dag_3_X.PitchTransferFcn_CSTATE +
    helicopter_dag_3_P.PitchTransferFcn_D * helicopter_dag_3_B.PitchCounttorad) *
    helicopter_dag_3_P.Gain_Gain_ae;
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' incorporates:
     *  Gain: '<S4>/Elevation_gain'
     */
    helicopter_dag_3_B.ElevationCounttorad = helicopter_dag_3_P.elevation_gain *
      rtb_DeadZonex * helicopter_dag_3_P.ElevationCounttorad_Gain;

    /* Gain: '<S9>/Gain' */
    helicopter_dag_3_B.Gain_e = helicopter_dag_3_P.Gain_Gain_lv *
      helicopter_dag_3_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter_dag_3_B.Sum = helicopter_dag_3_B.Gain_e +
      helicopter_dag_3_P.elavation_offsetdeg_Value;
  }

  /* Gain: '<S10>/Gain' incorporates:
   *  TransferFcn: '<S4>/Elevation: Transfer Fcn'
   */
  helicopter_dag_3_B.Gain_dg = (helicopter_dag_3_P.ElevationTransferFcn_C *
    helicopter_dag_3_X.ElevationTransferFcn_CSTATE +
    helicopter_dag_3_P.ElevationTransferFcn_D *
    helicopter_dag_3_B.ElevationCounttorad) * helicopter_dag_3_P.Gain_Gain_n;

  /* Gain: '<S2>/Gain1' */
  rtb_Gain1_idx_2 = helicopter_dag_3_P.Gain1_Gain * helicopter_dag_3_B.Gain_i;
  rtb_Gain1_idx_3 = helicopter_dag_3_P.Gain1_Gain * helicopter_dag_3_B.Gain_b;

  /* Sum: '<Root>/Subtract' incorporates:
   *  Constant: '<Root>/Constant'
   *  Gain: '<S2>/Gain1'
   */
  helicopter_dag_3_B.Subtract[0] = helicopter_dag_3_P.x0[0] -
    helicopter_dag_3_P.Gain1_Gain * helicopter_dag_3_B.Gain;
  helicopter_dag_3_B.Subtract[1] = helicopter_dag_3_P.x0[1] -
    helicopter_dag_3_P.Gain1_Gain * helicopter_dag_3_B.Gain_d;
  helicopter_dag_3_B.Subtract[2] = helicopter_dag_3_P.x0[2] - rtb_Gain1_idx_2;
  helicopter_dag_3_B.Subtract[3] = helicopter_dag_3_P.x0[3] - rtb_Gain1_idx_3;
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
    /* SignalConversion generated from: '<Root>/To Workspace' */
    helicopter_dag_3_B.TmpSignalConversionAtToWorkspac[0] =
      helicopter_dag_3_B.FromWorkspace[0];
    helicopter_dag_3_B.TmpSignalConversionAtToWorkspac[4] =
      helicopter_dag_3_B.Subtract[0];
    helicopter_dag_3_B.TmpSignalConversionAtToWorkspac[1] =
      helicopter_dag_3_B.FromWorkspace[1];
    helicopter_dag_3_B.TmpSignalConversionAtToWorkspac[5] =
      helicopter_dag_3_B.Subtract[1];
    helicopter_dag_3_B.TmpSignalConversionAtToWorkspac[2] =
      helicopter_dag_3_B.FromWorkspace[2];
    helicopter_dag_3_B.TmpSignalConversionAtToWorkspac[6] =
      helicopter_dag_3_B.Subtract[2];
    helicopter_dag_3_B.TmpSignalConversionAtToWorkspac[3] =
      helicopter_dag_3_B.FromWorkspace[3];
    helicopter_dag_3_B.TmpSignalConversionAtToWorkspac[7] =
      helicopter_dag_3_B.Subtract[3];
  }

  /* FromWorkspace: '<S5>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      helicopter_dag_3_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter_dag_3_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = helicopter_dag_3_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helicopter_dag_3_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicopter_dag_3_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter_dag_3_B.FromWorkspace1[0])[elIdx] =
                pDataValues[currTimeIndex];
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter_dag_3_B.FromWorkspace1[0])[elIdx] =
                pDataValues[currTimeIndex + 1];
              pDataValues += 141;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 4; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&helicopter_dag_3_B.FromWorkspace1[0])[elIdx] = (real_T)
              rtInterpolate(d1, d2, f1, f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

  /* FromWorkspace: '<S5>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helicopter_dag_3_DW.FromWorkspace_PWORK_e.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter_dag_3_DW.FromWorkspace_PWORK_e.TimePtr;
    int_T currTimeIndex = helicopter_dag_3_DW.FromWorkspace_IWORK_p.PrevIndex;
    real_T t = helicopter_dag_3_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicopter_dag_3_DW.FromWorkspace_IWORK_p.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Frontgain = pDataValues[currTimeIndex];
        } else {
          rtb_Frontgain = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Frontgain = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  /* Sum: '<S5>/Subtract' incorporates:
   *  Gain: '<S5>/Gain'
   *  Sum: '<S5>/Subtract1'
   */
  helicopter_dag_3_B.Subtract_a = rtb_Frontgain -
    ((((helicopter_dag_3_B.FromWorkspace1[0] - helicopter_dag_3_B.Subtract[0]) *
       helicopter_dag_3_P.K1[0] + (helicopter_dag_3_B.FromWorkspace1[1] -
        helicopter_dag_3_B.Subtract[1]) * helicopter_dag_3_P.K1[1]) +
      (helicopter_dag_3_B.FromWorkspace1[2] - helicopter_dag_3_B.Subtract[2]) *
      helicopter_dag_3_P.K1[2]) + (helicopter_dag_3_B.FromWorkspace1[3] -
      helicopter_dag_3_B.Subtract[3]) * helicopter_dag_3_P.K1[3]);
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S7>/Sum2'
   */
  rtb_Frontgain = helicopter_dag_3_B.Subtract_a - rtb_Gain1_idx_2;

  /* Gain: '<S7>/K_pp' */
  rtb_Gain1_idx_2 = helicopter_dag_3_P.K_pp * rtb_Frontgain;

  /* Gain: '<S1>/Front gain' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S7>/K_pd'
   *  Sum: '<Root>/Sum1'
   *  Sum: '<S7>/Sum3'
   */
  rtb_Frontgain = helicopter_dag_3_P.K_pd * rtb_Gain1_idx_3;
  rtb_Frontgain = (rtb_Gain1_idx_2 - rtb_Frontgain) + helicopter_dag_3_P.Vd_ff;

  /* Integrator: '<S3>/Integrator' */
  /* Limited  Integrator  */
  if (helicopter_dag_3_X.Integrator_CSTATE >=
      helicopter_dag_3_P.Integrator_UpperSat) {
    helicopter_dag_3_X.Integrator_CSTATE =
      helicopter_dag_3_P.Integrator_UpperSat;
  } else {
    if (helicopter_dag_3_X.Integrator_CSTATE <=
        helicopter_dag_3_P.Integrator_LowerSat) {
      helicopter_dag_3_X.Integrator_CSTATE =
        helicopter_dag_3_P.Integrator_LowerSat;
    }
  }

  /* Sum: '<S3>/Sum' incorporates:
   *  Constant: '<Root>/elevation_ref'
   *  Gain: '<S2>/Gain1'
   */
  rtb_Gain1_idx_3 = helicopter_dag_3_P.elevation_ref_Value -
    helicopter_dag_3_P.Gain1_Gain * helicopter_dag_3_B.Sum;

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Vs_bias'
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Integrator: '<S3>/Integrator'
   *  Sum: '<S3>/Sum1'
   */
  rtb_Backgain = ((helicopter_dag_3_P.K_ep * rtb_Gain1_idx_3 +
                   helicopter_dag_3_X.Integrator_CSTATE) -
                  helicopter_dag_3_P.Gain1_Gain * helicopter_dag_3_B.Gain_dg *
                  helicopter_dag_3_P.K_ed) + helicopter_dag_3_P.Vs_ff;

  /* Sum: '<S1>/Subtract' */
  rtb_Gain1_idx_2 = rtb_Backgain - rtb_Frontgain;

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Frontgain = (rtb_Frontgain + rtb_Backgain) *
    helicopter_dag_3_P.Frontgain_Gain;

  /* If: '<S3>/If' incorporates:
   *  Clock: '<S3>/Clock'
   *  Gain: '<S3>/K_ei'
   *  Inport: '<S8>/In1'
   */
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
    rtAction = (int8_T)!(helicopter_dag_3_M->Timing.t[0] >= 2.0);
    helicopter_dag_3_DW.If_ActiveSubsystem = rtAction;
  } else {
    rtAction = helicopter_dag_3_DW.If_ActiveSubsystem;
  }

  if (rtAction == 0) {
    /* Outputs for IfAction SubSystem: '<S3>/If Action Subsystem' incorporates:
     *  ActionPort: '<S8>/Action Port'
     */
    helicopter_dag_3_B.In1 = helicopter_dag_3_P.K_ei * rtb_Gain1_idx_3;
    if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
      srUpdateBC(helicopter_dag_3_DW.IfActionSubsystem_SubsysRanBC);
    }

    /* End of Outputs for SubSystem: '<S3>/If Action Subsystem' */
  }

  /* End of If: '<S3>/If' */
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  rtb_Gain1_idx_3 = helicopter_dag_3_M->Timing.t[0];
  if ((helicopter_dag_3_DW.TimeStampA >= rtb_Gain1_idx_3) &&
      (helicopter_dag_3_DW.TimeStampB >= rtb_Gain1_idx_3)) {
    rtb_Gain1_idx_3 = 0.0;
  } else {
    rtb_Backgain = helicopter_dag_3_DW.TimeStampA;
    lastU = &helicopter_dag_3_DW.LastUAtTimeA;
    if (helicopter_dag_3_DW.TimeStampA < helicopter_dag_3_DW.TimeStampB) {
      if (helicopter_dag_3_DW.TimeStampB < rtb_Gain1_idx_3) {
        rtb_Backgain = helicopter_dag_3_DW.TimeStampB;
        lastU = &helicopter_dag_3_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter_dag_3_DW.TimeStampA >= rtb_Gain1_idx_3) {
        rtb_Backgain = helicopter_dag_3_DW.TimeStampB;
        lastU = &helicopter_dag_3_DW.LastUAtTimeB;
      }
    }

    rtb_Gain1_idx_3 = (helicopter_dag_3_B.PitchCounttorad - *lastU) /
      (rtb_Gain1_idx_3 - rtb_Backgain);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S13>/Gain' */
  helicopter_dag_3_B.Gain_l = helicopter_dag_3_P.Gain_Gain_a1 * rtb_Gain1_idx_3;
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
  }

  /* Gain: '<S1>/Back gain' */
  rtb_Gain1_idx_2 *= helicopter_dag_3_P.Backgain_Gain;

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Gain1_idx_2 > helicopter_dag_3_P.BackmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Back motor: Saturation' */
    helicopter_dag_3_B.BackmotorSaturation =
      helicopter_dag_3_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Gain1_idx_2 < helicopter_dag_3_P.BackmotorSaturation_LowerSat)
  {
    /* Saturate: '<S4>/Back motor: Saturation' */
    helicopter_dag_3_B.BackmotorSaturation =
      helicopter_dag_3_P.BackmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Back motor: Saturation' */
    helicopter_dag_3_B.BackmotorSaturation = rtb_Gain1_idx_2;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
  }

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Frontgain > helicopter_dag_3_P.FrontmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Front motor: Saturation' */
    helicopter_dag_3_B.FrontmotorSaturation =
      helicopter_dag_3_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Frontgain < helicopter_dag_3_P.FrontmotorSaturation_LowerSat) {
    /* Saturate: '<S4>/Front motor: Saturation' */
    helicopter_dag_3_B.FrontmotorSaturation =
      helicopter_dag_3_P.FrontmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Front motor: Saturation' */
    helicopter_dag_3_B.FrontmotorSaturation = rtb_Frontgain;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helicopter_dag_3/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter_dag_3_DW.HILWriteAnalog_Buffer[0] =
        helicopter_dag_3_B.FrontmotorSaturation;
      helicopter_dag_3_DW.HILWriteAnalog_Buffer[1] =
        helicopter_dag_3_B.BackmotorSaturation;
      result = hil_write_analog(helicopter_dag_3_DW.HILInitialize_Card,
        helicopter_dag_3_P.HILWriteAnalog_channels, 2,
        &helicopter_dag_3_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
      }
    }

    /* S-Function (game_controller_block): '<S6>/Game Controller' */

    /* S-Function Block: helicopter_dag_3/Joystick/Game Controller (game_controller_block) */
    {
      if (helicopter_dag_3_P.GameController_Enabled) {
        t_game_controller_states state;
        t_boolean new_data;
        t_error result;
        result = game_controller_poll
          (helicopter_dag_3_DW.GameController_Controller, &state, &new_data);
        if (result < 0) {
          new_data = false;
        }

        rtb_DeadZonex = state.x;
        rtb_DeadZoney = state.y;
      } else {
        rtb_DeadZonex = 0;
        rtb_DeadZoney = 0;
      }
    }

    /* DeadZone: '<S6>/Dead Zone: x' */
    if (rtb_DeadZonex > helicopter_dag_3_P.DeadZonex_End) {
      /* DeadZone: '<S6>/Dead Zone: x' */
      rtb_DeadZonex -= helicopter_dag_3_P.DeadZonex_End;
    } else if (rtb_DeadZonex >= helicopter_dag_3_P.DeadZonex_Start) {
      /* DeadZone: '<S6>/Dead Zone: x' */
      rtb_DeadZonex = 0.0;
    } else {
      /* DeadZone: '<S6>/Dead Zone: x' */
      rtb_DeadZonex -= helicopter_dag_3_P.DeadZonex_Start;
    }

    /* End of DeadZone: '<S6>/Dead Zone: x' */

    /* Gain: '<S6>/Joystick_gain_x' incorporates:
     *  Gain: '<S6>/Gain: x'
     */
    helicopter_dag_3_B.Joystick_gain_x = helicopter_dag_3_P.Gainx_Gain *
      rtb_DeadZonex * helicopter_dag_3_P.Joystick_gain_x_Gain;

    /* DeadZone: '<S6>/Dead Zone: y' */
    if (rtb_DeadZoney > helicopter_dag_3_P.DeadZoney_End) {
      /* DeadZone: '<S6>/Dead Zone: y' */
      rtb_DeadZoney -= helicopter_dag_3_P.DeadZoney_End;
    } else if (rtb_DeadZoney >= helicopter_dag_3_P.DeadZoney_Start) {
      /* DeadZone: '<S6>/Dead Zone: y' */
      rtb_DeadZoney = 0.0;
    } else {
      /* DeadZone: '<S6>/Dead Zone: y' */
      rtb_DeadZoney -= helicopter_dag_3_P.DeadZoney_Start;
    }

    /* End of DeadZone: '<S6>/Dead Zone: y' */

    /* Gain: '<S6>/Joystick_gain_y' incorporates:
     *  Gain: '<S6>/Gain: y'
     */
    helicopter_dag_3_B.Joystick_gain_y = helicopter_dag_3_P.Gainy_Gain *
      rtb_DeadZoney * helicopter_dag_3_P.Joystick_gain_y_Gain;
  }
}

/* Model update function */
void helicopter_dag_3_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (helicopter_dag_3_DW.TimeStampA == (rtInf)) {
    helicopter_dag_3_DW.TimeStampA = helicopter_dag_3_M->Timing.t[0];
    lastU = &helicopter_dag_3_DW.LastUAtTimeA;
  } else if (helicopter_dag_3_DW.TimeStampB == (rtInf)) {
    helicopter_dag_3_DW.TimeStampB = helicopter_dag_3_M->Timing.t[0];
    lastU = &helicopter_dag_3_DW.LastUAtTimeB;
  } else if (helicopter_dag_3_DW.TimeStampA < helicopter_dag_3_DW.TimeStampB) {
    helicopter_dag_3_DW.TimeStampA = helicopter_dag_3_M->Timing.t[0];
    lastU = &helicopter_dag_3_DW.LastUAtTimeA;
  } else {
    helicopter_dag_3_DW.TimeStampB = helicopter_dag_3_M->Timing.t[0];
    lastU = &helicopter_dag_3_DW.LastUAtTimeB;
  }

  *lastU = helicopter_dag_3_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter_dag_3_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helicopter_dag_3_M->Timing.clockTick0)) {
    ++helicopter_dag_3_M->Timing.clockTickH0;
  }

  helicopter_dag_3_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicopter_dag_3_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.02s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helicopter_dag_3_M->Timing.clockTick1)) {
      ++helicopter_dag_3_M->Timing.clockTickH1;
    }

    helicopter_dag_3_M->Timing.t[1] = helicopter_dag_3_M->Timing.clockTick1 *
      helicopter_dag_3_M->Timing.stepSize1 +
      helicopter_dag_3_M->Timing.clockTickH1 *
      helicopter_dag_3_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter_dag_3_derivatives(void)
{
  XDot_helicopter_dag_3_T *_rtXdot;
  boolean_T lsat;
  boolean_T usat;
  _rtXdot = ((XDot_helicopter_dag_3_T *) helicopter_dag_3_M->derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter_dag_3_P.TravelTransferFcn_A *
    helicopter_dag_3_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter_dag_3_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_dag_3_P.PitchTransferFcn_A *
    helicopter_dag_3_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_dag_3_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter_dag_3_P.ElevationTransferFcn_A *
    helicopter_dag_3_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter_dag_3_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  lsat = (helicopter_dag_3_X.Integrator_CSTATE <=
          helicopter_dag_3_P.Integrator_LowerSat);
  usat = (helicopter_dag_3_X.Integrator_CSTATE >=
          helicopter_dag_3_P.Integrator_UpperSat);
  if (((!lsat) && (!usat)) || (lsat && (helicopter_dag_3_B.In1 > 0.0)) || (usat &&
       (helicopter_dag_3_B.In1 < 0.0))) {
    _rtXdot->Integrator_CSTATE = helicopter_dag_3_B.In1;
  } else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S3>/Integrator' */
}

/* Model initialize function */
void helicopter_dag_3_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter_dag_3/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter_dag_3_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options
      (helicopter_dag_3_DW.HILInitialize_Card, "update_rate=normal;decimation=1",
       32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter_dag_3_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
      return;
    }

    if ((helicopter_dag_3_P.HILInitialize_AIPStart && !is_switching) ||
        (helicopter_dag_3_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter_dag_3_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = (helicopter_dag_3_P.HILInitialize_AILow);
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter_dag_3_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helicopter_dag_3_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges
        (helicopter_dag_3_DW.HILInitialize_Card,
         helicopter_dag_3_P.HILInitialize_AIChannels, 8U,
         &helicopter_dag_3_DW.HILInitialize_AIMinimums[0],
         &helicopter_dag_3_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_dag_3_P.HILInitialize_AOPStart && !is_switching) ||
        (helicopter_dag_3_P.HILInitialize_AOPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter_dag_3_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = (helicopter_dag_3_P.HILInitialize_AOLow);
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter_dag_3_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helicopter_dag_3_P.HILInitialize_AOHigh;
        }
      }

      result = hil_set_analog_output_ranges
        (helicopter_dag_3_DW.HILInitialize_Card,
         helicopter_dag_3_P.HILInitialize_AOChannels, 8U,
         &helicopter_dag_3_DW.HILInitialize_AOMinimums[0],
         &helicopter_dag_3_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_dag_3_P.HILInitialize_AOStart && !is_switching) ||
        (helicopter_dag_3_P.HILInitialize_AOEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_dag_3_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_dag_3_P.HILInitialize_AOInitial;
        }
      }

      result = hil_write_analog(helicopter_dag_3_DW.HILInitialize_Card,
        helicopter_dag_3_P.HILInitialize_AOChannels, 8U,
        &helicopter_dag_3_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_dag_3_P.HILInitialize_AOReset) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_dag_3_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_dag_3_P.HILInitialize_AOWatchdog;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter_dag_3_DW.HILInitialize_Card,
         helicopter_dag_3_P.HILInitialize_AOChannels, 8U,
         &helicopter_dag_3_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_dag_3_P.HILInitialize_EIPStart && !is_switching) ||
        (helicopter_dag_3_P.HILInitialize_EIPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter_dag_3_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter_dag_3_P.HILInitialize_EIQuadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicopter_dag_3_DW.HILInitialize_Card,
         helicopter_dag_3_P.HILInitialize_EIChannels, 8U,
         (t_encoder_quadrature_mode *)
         &helicopter_dag_3_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_dag_3_P.HILInitialize_EIStart && !is_switching) ||
        (helicopter_dag_3_P.HILInitialize_EIEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter_dag_3_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = helicopter_dag_3_P.HILInitialize_EIInitial;
        }
      }

      result = hil_set_encoder_counts(helicopter_dag_3_DW.HILInitialize_Card,
        helicopter_dag_3_P.HILInitialize_EIChannels, 8U,
        &helicopter_dag_3_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_dag_3_P.HILInitialize_POPStart && !is_switching) ||
        (helicopter_dag_3_P.HILInitialize_POPEnter && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter_dag_3_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter_dag_3_P.HILInitialize_POModes;
        }
      }

      result = hil_set_pwm_mode(helicopter_dag_3_DW.HILInitialize_Card,
        helicopter_dag_3_P.HILInitialize_POChannels, 8U, (t_pwm_mode *)
        &helicopter_dag_3_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_POChannels =
          helicopter_dag_3_P.HILInitialize_POChannels;
        int32_T *dw_POModeValues =
          &helicopter_dag_3_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE ||
              dw_POModeValues[i1] == PWM_RAW_MODE) {
            helicopter_dag_3_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = (p_HILInitialize_POChannels[i1]);
            helicopter_dag_3_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = helicopter_dag_3_P.HILInitialize_POFrequency;
            num_duty_cycle_modes++;
          } else {
            helicopter_dag_3_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = (p_HILInitialize_POChannels[i1]);
            helicopter_dag_3_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] =
              helicopter_dag_3_P.HILInitialize_POFrequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter_dag_3_DW.HILInitialize_Card,
          &helicopter_dag_3_DW.HILInitialize_POSortedChans[0],
          num_duty_cycle_modes,
          &helicopter_dag_3_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter_dag_3_DW.HILInitialize_Card,
          &helicopter_dag_3_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter_dag_3_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter_dag_3_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter_dag_3_P.HILInitialize_POConfiguration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicopter_dag_3_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicopter_dag_3_P.HILInitialize_POAlignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter_dag_3_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicopter_dag_3_P.HILInitialize_POPolarity;
        }
      }

      result = hil_set_pwm_configuration(helicopter_dag_3_DW.HILInitialize_Card,
        helicopter_dag_3_P.HILInitialize_POChannels, 8U,
        (t_pwm_configuration *) &helicopter_dag_3_DW.HILInitialize_POModeValues
        [0],
        (t_pwm_alignment *) &helicopter_dag_3_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicopter_dag_3_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &helicopter_dag_3_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = helicopter_dag_3_P.HILInitialize_POLeading;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter_dag_3_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_dag_3_P.HILInitialize_POTrailing;
        }
      }

      result = hil_set_pwm_deadband(helicopter_dag_3_DW.HILInitialize_Card,
        helicopter_dag_3_P.HILInitialize_POChannels, 8U,
        &helicopter_dag_3_DW.HILInitialize_POSortedFreqs[0],
        &helicopter_dag_3_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_dag_3_P.HILInitialize_POStart && !is_switching) ||
        (helicopter_dag_3_P.HILInitialize_POEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_dag_3_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_dag_3_P.HILInitialize_POInitial;
        }
      }

      result = hil_write_pwm(helicopter_dag_3_DW.HILInitialize_Card,
        helicopter_dag_3_P.HILInitialize_POChannels, 8U,
        &helicopter_dag_3_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_dag_3_P.HILInitialize_POReset) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_dag_3_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_dag_3_P.HILInitialize_POWatchdog;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter_dag_3_DW.HILInitialize_Card,
         helicopter_dag_3_P.HILInitialize_POChannels, 8U,
         &helicopter_dag_3_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter_dag_3/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (helicopter_dag_3_DW.HILInitialize_Card,
       helicopter_dag_3_P.HILReadEncoderTimebase_SamplesI,
       helicopter_dag_3_P.HILReadEncoderTimebase_Channels, 3,
       &helicopter_dag_3_DW.HILReadEncoderTimebase_Task);
    if (result >= 0) {
      result = hil_task_set_buffer_overflow_mode
        (helicopter_dag_3_DW.HILReadEncoderTimebase_Task,
         (t_buffer_overflow_mode)
         (helicopter_dag_3_P.HILReadEncoderTimebase_Overflow - 1));
    }

    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
    }
  }

  /* Start for FromWorkspace: '<Root>/From Workspace' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1378421413625257, 3.1262155534579978,
      3.1033093000299643, 3.0666274151911783, 3.0144539223941584,
      2.9456562771175667, 2.8595077632935446, 2.7555515879651526,
      2.633505110490284, 2.4931956060320961, 2.334518576064299,
      2.1574113214711184, 1.9626074078533053, 1.7564316060085343,
      1.5468504786963275, 1.3407140017858634, 1.1431472634984829,
      0.95771574317723407, 0.78679112891251446, 0.63190372010248919,
      0.49401941043779485, 0.37373776741018649, 0.27142613241917724,
      0.18730721183494006, 0.12060531290453495, 0.069687620981889756,
      0.032423160561046827, 0.0064951272092745194, -0.010374343815386126,
      -0.020279599214938396, -0.025050584982836358, -0.026219158735923551,
      -0.025018303281833194, -0.022402620880901388, -0.019081092484360112,
      -0.015555420158210679, -0.012159253593134033, -0.0090952137019660977,
      -0.0064678838098550188, -0.0043118781455844732, -0.0026147646291379391,
      -0.0013350637243603164, -0.00041581512383159438, 0.00020565753593144271,
      0.00059110968550767323, 0.000797212021613826, 0.00087326458104617592,
      0.00086027374111505577, 0.00079096539499473042, 0.00069039448079168076,
      0.00057689096058666429, 0.00046315206552779238, 0.00035734883159368227,
      0.00026416162815334239, 0.00018569533921450638, 0.00012225138129625049,
      7.2952330495698411E-5, 3.6227111707429084E-5, 1.0171911483564643E-5,
      -7.19451186936615E-6, -1.776155061641052E-5, -2.322477144568408E-5,
      -2.50343986275888E-5, -2.4379837501495205E-5, -2.2198409578408559E-5,
      -1.9198980842800539E-5, -1.5893442986654058E-5, -1.2630972733833013E-5,
      -9.6316181137728373E-6, -7.0170485935265626E-6, -4.8372864436166657E-6,
      -3.0929497086123142E-6, -1.7530276070961006E-6, -7.6852187299918131E-7,
      -8.2463923703578443E-8, 3.6310585756364069E-7, 6.226035502259896E-7,
      7.4425679063593389E-7, 7.6885871231037426E-7, 7.2950624230653111E-7,
      6.5198686983039307E-7, 5.5555242755263239E-7, 4.5388483372047351E-7,
      3.5611525313438853E-7, 2.6780425257410815E-7, 1.918266792641321E-7,
      1.2913216679472209E-7, 7.9371645807404458E-8, 4.1393372434475351E-8,
      1.3620111282669631E-8, -5.6765744340085384E-9, -1.8187426373784671E-8,
      -2.544485342868394E-8, -2.875501341395178E-8, -2.9169383662956152E-8,
      -2.7485240964150037E-8, -2.4267230479563931E-8, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048909068408,
      -0.046506351618112091, -0.091625013712135384, -0.14672753935514368,
      -0.20869397118807922, -0.27519058110636674, -0.34459405529608839,
      -0.41582470131356875, -0.48818590989947475, -0.561238017832752,
      -0.63470811987118869, -0.70842901837272132, -0.77921565447125174,
      -0.82470320737908487, -0.8383245092488274, -0.8245459076418562,
      -0.79026695314952233, -0.74172608128499551, -0.68369845705887855,
      -0.61954963524010087, -0.55153723865877735, -0.48112657211043336,
      -0.409246539964037, -0.33647568233694863, -0.26680759572162049,
      -0.20367076769058079, -0.14905784168337169, -0.10371213340708924,
      -0.067477884098642582, -0.039621021598209079, -0.019083943071591846,
      -0.0046742950123487515, 0.0048034218163614155, 0.010462729603727221,
      0.0132861135861651, 0.014102689304597732, 0.013584666260306582,
      0.012256159564671744, 0.010509319568444312, 0.0086240226570821858,
      0.0067884540657861354, 0.0051188036191104911, 0.0036769944021148877,
      0.0024858906390521486, 0.0015418085983049221, 0.00082440934442461124,
      0.00030421023772939982, -5.1963359724480886E-5, -0.00027723338448130123,
      -0.0004022836568121983, -0.00045401408082006563, -0.00045495558023548785,
      -0.0004232129357364406, -0.00037274881376135957, -0.00031386515575534393,
      -0.0002537758316730235, -0.0001971962032022083, -0.00014690087515307731,
      -0.00010422080089545776, -6.9465693411723182E-5, -4.226815498817748E-5,
      -2.1852883317094237E-5, -7.2385087276188866E-6, 2.6182445043743876E-6,
      8.7257116923465777E-6, 1.1997714942432086E-5, 1.322215142458592E-5,
      1.3049881011284181E-5, 1.1997418480240705E-5, 1.0458278080985097E-5,
      8.7190485996395858E-6, 6.9773469400174062E-6, 5.3596884060648541E-6,
      3.9380229363876772E-6, 2.7442317971824115E-6, 1.7822791250688765E-6,
      1.0379907706493961E-6, 4.8661296163977715E-7, 9.8407686697761432E-8,
      -1.5740988001537237E-7, -3.1007748990455181E-7, -3.85737769111043E-7,
      -4.0667037532863545E-7, -3.9107832234433978E-7, -3.532440022411216E-7,
      -3.0391029323990427E-7, -2.5077804987764E-7, -1.9904208394927056E-7,
      -1.519130934917164E-7, -1.1109304460722289E-7, -7.7186742866712678E-8,
      -5.0043407759104529E-8, -2.9029708219597069E-8, -1.3240639941071359E-8,
      -1.6574809960174952E-9, 6.7365707952244546E-9, 1.287204193834443E-8,
      1.7528792604517883E-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875205865547, 0.2226603793231765, 0.31888147181640647,
      0.3894436063114417, 0.43795507377677845, 0.46997264230390068,
      0.49051724877547087, 0.50343100141474351, 0.51142138586029351,
      0.51630439857701849, 0.519258621270637, 0.52103115488680818,
      0.50029290888529865, 0.321488651200451, 0.096270159324037929,
      -0.097381893790350116, -0.2422705584235621, -0.34306834345299281,
      -0.41011708593382, -0.45337937269132711, -0.48068564352413468,
      -0.497635699696077, -0.50802061455828917, -0.51431662882889484,
      -0.49238742833803084, -0.44622698710350855, -0.385983303075076,
      -0.32048614751840826, -0.2560898354113797, -0.19688166497078385,
      -0.14514822743922479, -0.10184188910372288, -0.0669848827853115,
      -0.039997826020118532, -0.019954599672004081, -0.0057712453086642412,
      0.0036611890320270568, 0.00938937793720851, 0.012345997934600894,
      0.013324558530863784, 0.012973097757036256, 0.011800451678849244,
      0.010190156885323987, 0.0084182664871576, 0.006672411297160119,
      0.0050703039350017409, 0.0036765686100103112, 0.0025172989558410741,
      0.0015921224991307659, 0.00088380756523365189, 0.00036561087983799645,
      6.6541582877110983E-6, -0.0002243448880666632, -0.00035666114070065369,
      -0.00041616720575099464, -0.0004246883931064982, -0.00039988320495853991,
      -0.00035546816969822537, -0.00030164646434205356, -0.00024563582591297095,
      -0.00019222181420652795, -0.00014428734310512059, -0.00010328881802967871,
      -6.9663767320982117E-5, -4.3165245501186611E-5, -2.3125269317691988E-5,
      -8.6538494151744416E-6, 1.2175414870752022E-6, 7.4384032094920727E-6,
      1.0878056508456524E-5, 1.229220972209788E-5, 1.2309682133948563E-5,
      1.14329926965695E-5, 1.0047788572609129E-5, 8.43725983534771E-6,
      6.798714094480296E-6, 5.2603458279287807E-6, 3.8969277701417226E-6,
      2.7436866187890985E-6, 1.808020858917736E-6, 1.0789963593715868E-6,
      5.3473795691427739E-7, 1.4794366609560683E-7, -1.101986755180917E-7,
      -2.6739852465595959E-7, -3.4867181364361954E-7, -3.7551840370575462E-7,
      -3.6565004801758505E-7, -3.3308970492917922E-7, -2.8850051536633003E-7,
      -2.3963679118832459E-7, -1.9183872590478046E-7, -1.48516802744858E-7,
      -1.1159110435698949E-7, -8.1865343504183841E-8, -5.932595215085712E-8,
      -4.3363166724041946E-8, -3.291213523937131E-8, -2.6497511518286387E-8, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500823462206,
      0.46652650905808429, 0.38488436997291947, 0.282248537980141,
      0.19404586986134692, 0.12807027410848898, 0.082178425886280229,
      0.051655010557090514, 0.031961537782199671, 0.019532050866899832,
      0.011816890774474108, 0.0070901344646842534, -0.082952984006037914,
      -0.71521703073939069, -0.90087396750565218, -0.77460821245755218,
      -0.57955465853284782, -0.40319114011772295, -0.26819496992330877,
      -0.17304914703002844, -0.10922508333123025, -0.067800224687769459,
      -0.0415396594488484, -0.025184057082422755, 0.087716801963455965,
      0.18464176493808909, 0.24097473611373016, 0.26198862222667119,
      0.25758524842811414, 0.23683268176238331, 0.20693375012623633,
      0.1732253533420077, 0.13942802527364562, 0.10794822706077176,
      0.080172905392457749, 0.056733417453359435, 0.0377297373627655,
      0.022912755620725721, 0.01182647998956976, 0.00391424238505147,
      -0.0014058430953100788, -0.0046905843127481844, -0.0064411791741008553,
      -0.007087561592665317, -0.0069834207599896855, -0.0064084294486337866,
      -0.0055749412999657466, -0.0046370786166767351, -0.0037007058268414216,
      -0.002833259735588063, -0.0020727867415825402, -0.0014358268862016449,
      -0.00092399618541722124, -0.00052926501053607918, -0.00023802426020138118,
      -3.4084749422378092E-5, 9.9220752591737562E-5, 0.00017766014104137891,
      0.00021528682142488075, 0.00022404255371599979, 0.0002136560468258419,
      0.00019173788440595785, 0.00016399410030173341, 0.00013450020283461834,
      0.00010599408727900231, 8.01599047343207E-5, 5.7885679609559545E-5,
      3.9485563609347194E-5, 2.4883446889708624E-5, 1.3758613195885691E-5,
      5.6566128546446222E-6, 6.9889647437555716E-8, -3.5067577492915047E-6,
      -5.5408164961564335E-6, -6.4421149489851418E-6, -6.5541829636083467E-6,
      -6.1534730660817039E-6, -5.4536722316848E-6, -4.61296460560862E-6,
      -3.7426630393699866E-6, -2.9160979980771008E-6, -2.1770336094460068E-6,
      -1.5471771634872629E-6, -1.0325693660781116E-6, -6.2879939682010233E-7,
      -3.2509315591580272E-7, -1.073863603259355E-7, 3.9473422770938195E-8,
      1.3024137263752579E-7, 1.783567579302375E-7, 1.9545489689867386E-7,
      1.9119226140673845E-7, 1.7328769237973729E-7, 1.4770279392306797E-7,
      1.1890304314760869E-7, 9.0157565566825517E-8, 6.3851141971350387E-8,
      4.1804125937105562E-8, 2.56584945530201E-8, 1.7646687417935339E-8, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0 } ;

    helicopter_dag_3_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter_dag_3_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helicopter_dag_3_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<S5>/From Workspace1' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1378421413625257, 3.1262155534579978,
      3.1033093000299643, 3.0666274151911783, 3.0144539223941584,
      2.9456562771175667, 2.8595077632935446, 2.7555515879651526,
      2.633505110490284, 2.4931956060320961, 2.334518576064299,
      2.1574113214711184, 1.9626074078533053, 1.7564316060085343,
      1.5468504786963275, 1.3407140017858634, 1.1431472634984829,
      0.95771574317723407, 0.78679112891251446, 0.63190372010248919,
      0.49401941043779485, 0.37373776741018649, 0.27142613241917724,
      0.18730721183494006, 0.12060531290453495, 0.069687620981889756,
      0.032423160561046827, 0.0064951272092745194, -0.010374343815386126,
      -0.020279599214938396, -0.025050584982836358, -0.026219158735923551,
      -0.025018303281833194, -0.022402620880901388, -0.019081092484360112,
      -0.015555420158210679, -0.012159253593134033, -0.0090952137019660977,
      -0.0064678838098550188, -0.0043118781455844732, -0.0026147646291379391,
      -0.0013350637243603164, -0.00041581512383159438, 0.00020565753593144271,
      0.00059110968550767323, 0.000797212021613826, 0.00087326458104617592,
      0.00086027374111505577, 0.00079096539499473042, 0.00069039448079168076,
      0.00057689096058666429, 0.00046315206552779238, 0.00035734883159368227,
      0.00026416162815334239, 0.00018569533921450638, 0.00012225138129625049,
      7.2952330495698411E-5, 3.6227111707429084E-5, 1.0171911483564643E-5,
      -7.19451186936615E-6, -1.776155061641052E-5, -2.322477144568408E-5,
      -2.50343986275888E-5, -2.4379837501495205E-5, -2.2198409578408559E-5,
      -1.9198980842800539E-5, -1.5893442986654058E-5, -1.2630972733833013E-5,
      -9.6316181137728373E-6, -7.0170485935265626E-6, -4.8372864436166657E-6,
      -3.0929497086123142E-6, -1.7530276070961006E-6, -7.6852187299918131E-7,
      -8.2463923703578443E-8, 3.6310585756364069E-7, 6.226035502259896E-7,
      7.4425679063593389E-7, 7.6885871231037426E-7, 7.2950624230653111E-7,
      6.5198686983039307E-7, 5.5555242755263239E-7, 4.5388483372047351E-7,
      3.5611525313438853E-7, 2.6780425257410815E-7, 1.918266792641321E-7,
      1.2913216679472209E-7, 7.9371645807404458E-8, 4.1393372434475351E-8,
      1.3620111282669631E-8, -5.6765744340085384E-9, -1.8187426373784671E-8,
      -2.544485342868394E-8, -2.875501341395178E-8, -2.9169383662956152E-8,
      -2.7485240964150037E-8, -2.4267230479563931E-8, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048909068408,
      -0.046506351618112091, -0.091625013712135384, -0.14672753935514368,
      -0.20869397118807922, -0.27519058110636674, -0.34459405529608839,
      -0.41582470131356875, -0.48818590989947475, -0.561238017832752,
      -0.63470811987118869, -0.70842901837272132, -0.77921565447125174,
      -0.82470320737908487, -0.8383245092488274, -0.8245459076418562,
      -0.79026695314952233, -0.74172608128499551, -0.68369845705887855,
      -0.61954963524010087, -0.55153723865877735, -0.48112657211043336,
      -0.409246539964037, -0.33647568233694863, -0.26680759572162049,
      -0.20367076769058079, -0.14905784168337169, -0.10371213340708924,
      -0.067477884098642582, -0.039621021598209079, -0.019083943071591846,
      -0.0046742950123487515, 0.0048034218163614155, 0.010462729603727221,
      0.0132861135861651, 0.014102689304597732, 0.013584666260306582,
      0.012256159564671744, 0.010509319568444312, 0.0086240226570821858,
      0.0067884540657861354, 0.0051188036191104911, 0.0036769944021148877,
      0.0024858906390521486, 0.0015418085983049221, 0.00082440934442461124,
      0.00030421023772939982, -5.1963359724480886E-5, -0.00027723338448130123,
      -0.0004022836568121983, -0.00045401408082006563, -0.00045495558023548785,
      -0.0004232129357364406, -0.00037274881376135957, -0.00031386515575534393,
      -0.0002537758316730235, -0.0001971962032022083, -0.00014690087515307731,
      -0.00010422080089545776, -6.9465693411723182E-5, -4.226815498817748E-5,
      -2.1852883317094237E-5, -7.2385087276188866E-6, 2.6182445043743876E-6,
      8.7257116923465777E-6, 1.1997714942432086E-5, 1.322215142458592E-5,
      1.3049881011284181E-5, 1.1997418480240705E-5, 1.0458278080985097E-5,
      8.7190485996395858E-6, 6.9773469400174062E-6, 5.3596884060648541E-6,
      3.9380229363876772E-6, 2.7442317971824115E-6, 1.7822791250688765E-6,
      1.0379907706493961E-6, 4.8661296163977715E-7, 9.8407686697761432E-8,
      -1.5740988001537237E-7, -3.1007748990455181E-7, -3.85737769111043E-7,
      -4.0667037532863545E-7, -3.9107832234433978E-7, -3.532440022411216E-7,
      -3.0391029323990427E-7, -2.5077804987764E-7, -1.9904208394927056E-7,
      -1.519130934917164E-7, -1.1109304460722289E-7, -7.7186742866712678E-8,
      -5.0043407759104529E-8, -2.9029708219597069E-8, -1.3240639941071359E-8,
      -1.6574809960174952E-9, 6.7365707952244546E-9, 1.287204193834443E-8,
      1.7528792604517883E-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875205865547, 0.2226603793231765, 0.31888147181640647,
      0.3894436063114417, 0.43795507377677845, 0.46997264230390068,
      0.49051724877547087, 0.50343100141474351, 0.51142138586029351,
      0.51630439857701849, 0.519258621270637, 0.52103115488680818,
      0.50029290888529865, 0.321488651200451, 0.096270159324037929,
      -0.097381893790350116, -0.2422705584235621, -0.34306834345299281,
      -0.41011708593382, -0.45337937269132711, -0.48068564352413468,
      -0.497635699696077, -0.50802061455828917, -0.51431662882889484,
      -0.49238742833803084, -0.44622698710350855, -0.385983303075076,
      -0.32048614751840826, -0.2560898354113797, -0.19688166497078385,
      -0.14514822743922479, -0.10184188910372288, -0.0669848827853115,
      -0.039997826020118532, -0.019954599672004081, -0.0057712453086642412,
      0.0036611890320270568, 0.00938937793720851, 0.012345997934600894,
      0.013324558530863784, 0.012973097757036256, 0.011800451678849244,
      0.010190156885323987, 0.0084182664871576, 0.006672411297160119,
      0.0050703039350017409, 0.0036765686100103112, 0.0025172989558410741,
      0.0015921224991307659, 0.00088380756523365189, 0.00036561087983799645,
      6.6541582877110983E-6, -0.0002243448880666632, -0.00035666114070065369,
      -0.00041616720575099464, -0.0004246883931064982, -0.00039988320495853991,
      -0.00035546816969822537, -0.00030164646434205356, -0.00024563582591297095,
      -0.00019222181420652795, -0.00014428734310512059, -0.00010328881802967871,
      -6.9663767320982117E-5, -4.3165245501186611E-5, -2.3125269317691988E-5,
      -8.6538494151744416E-6, 1.2175414870752022E-6, 7.4384032094920727E-6,
      1.0878056508456524E-5, 1.229220972209788E-5, 1.2309682133948563E-5,
      1.14329926965695E-5, 1.0047788572609129E-5, 8.43725983534771E-6,
      6.798714094480296E-6, 5.2603458279287807E-6, 3.8969277701417226E-6,
      2.7436866187890985E-6, 1.808020858917736E-6, 1.0789963593715868E-6,
      5.3473795691427739E-7, 1.4794366609560683E-7, -1.101986755180917E-7,
      -2.6739852465595959E-7, -3.4867181364361954E-7, -3.7551840370575462E-7,
      -3.6565004801758505E-7, -3.3308970492917922E-7, -2.8850051536633003E-7,
      -2.3963679118832459E-7, -1.9183872590478046E-7, -1.48516802744858E-7,
      -1.1159110435698949E-7, -8.1865343504183841E-8, -5.932595215085712E-8,
      -4.3363166724041946E-8, -3.291213523937131E-8, -2.6497511518286387E-8, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500823462206,
      0.46652650905808429, 0.38488436997291947, 0.282248537980141,
      0.19404586986134692, 0.12807027410848898, 0.082178425886280229,
      0.051655010557090514, 0.031961537782199671, 0.019532050866899832,
      0.011816890774474108, 0.0070901344646842534, -0.082952984006037914,
      -0.71521703073939069, -0.90087396750565218, -0.77460821245755218,
      -0.57955465853284782, -0.40319114011772295, -0.26819496992330877,
      -0.17304914703002844, -0.10922508333123025, -0.067800224687769459,
      -0.0415396594488484, -0.025184057082422755, 0.087716801963455965,
      0.18464176493808909, 0.24097473611373016, 0.26198862222667119,
      0.25758524842811414, 0.23683268176238331, 0.20693375012623633,
      0.1732253533420077, 0.13942802527364562, 0.10794822706077176,
      0.080172905392457749, 0.056733417453359435, 0.0377297373627655,
      0.022912755620725721, 0.01182647998956976, 0.00391424238505147,
      -0.0014058430953100788, -0.0046905843127481844, -0.0064411791741008553,
      -0.007087561592665317, -0.0069834207599896855, -0.0064084294486337866,
      -0.0055749412999657466, -0.0046370786166767351, -0.0037007058268414216,
      -0.002833259735588063, -0.0020727867415825402, -0.0014358268862016449,
      -0.00092399618541722124, -0.00052926501053607918, -0.00023802426020138118,
      -3.4084749422378092E-5, 9.9220752591737562E-5, 0.00017766014104137891,
      0.00021528682142488075, 0.00022404255371599979, 0.0002136560468258419,
      0.00019173788440595785, 0.00016399410030173341, 0.00013450020283461834,
      0.00010599408727900231, 8.01599047343207E-5, 5.7885679609559545E-5,
      3.9485563609347194E-5, 2.4883446889708624E-5, 1.3758613195885691E-5,
      5.6566128546446222E-6, 6.9889647437555716E-8, -3.5067577492915047E-6,
      -5.5408164961564335E-6, -6.4421149489851418E-6, -6.5541829636083467E-6,
      -6.1534730660817039E-6, -5.4536722316848E-6, -4.61296460560862E-6,
      -3.7426630393699866E-6, -2.9160979980771008E-6, -2.1770336094460068E-6,
      -1.5471771634872629E-6, -1.0325693660781116E-6, -6.2879939682010233E-7,
      -3.2509315591580272E-7, -1.073863603259355E-7, 3.9473422770938195E-8,
      1.3024137263752579E-7, 1.783567579302375E-7, 1.9545489689867386E-7,
      1.9119226140673845E-7, 1.7328769237973729E-7, 1.4770279392306797E-7,
      1.1890304314760869E-7, 9.0157565566825517E-8, 6.3851141971350387E-8,
      4.1804125937105562E-8, 2.56584945530201E-8, 1.7646687417935339E-8, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0 } ;

    helicopter_dag_3_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter_dag_3_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    helicopter_dag_3_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<S5>/From Workspace' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829859, 0.52359877559829859, 0.41597220466260454,
      -0.35171172454379279, -0.52359877559829771, -0.52359877559829837,
      -0.52359877559829848, -0.52359877559829848, -0.52359877559829848,
      -0.52359877559829837, -0.52359877559829826, -0.523598775598298,
      -0.52359877559829693, -0.52359877517451892, -0.39661912360557527,
      -0.297193067419335, -0.21768303374553549, -0.15253421103523923,
      -0.10032109788316607, -0.05990076860535648, -0.029854672510147673,
      -0.0085705806135195628, 0.0055992910045342681, 0.014214190567150231,
      0.018667071148491643, 0.02014554054051565, 0.019621197386550993,
      0.017858115042507561, 0.015433046349924817, 0.012761716685473323,
      0.010127149374335231, 0.0077072745762747363, 0.0056001091856460183,
      0.003845582943966197, 0.002443657378258357, 0.001368775903657804,
      0.0005809318464164992, 3.3779877371964417E-5, -0.00031972517415479462,
      -0.00052367530700314635, -0.00061708832558993176, -0.00063291862247416564,
      -0.00059786257299787593, -0.00053268583182874973, -0.0004528606401009716,
      -0.00036935536712312533, -0.0002894644557089876, -0.00021760436127671312,
      -0.00015603035641453822, -0.00010545104429882102, -6.553314155910428E-5,
      -3.5299678106448518E-5, -1.3431305735944754E-5, 1.5161048014800116E-6,
      1.0963116532702877E-5, 1.6213388242536375E-5, 1.8402271950823668E-5,
      1.847595987791717E-5, 1.7191694447471484E-5, 1.5131453224070945E-5,
      1.2723281647919471E-5, 1.0265992757774178E-5, 7.9542508157892655E-6,
      5.9021009966020443E-6, 4.1638157838308487E-6, 2.7515275000666861E-6,
      1.6495378974967778E-6, 8.2547344626426877E-7, 2.3862066445001773E-7,
      -1.5414208032193955E-7, -3.9437720111834551E-7, -5.196628567460948E-7,
      -5.6232932288402537E-7, -5.4902951196478256E-7, -5.0087418612676515E-7,
      -4.3391869553044415E-7, -3.5983932067384927E-7, -2.8668161755618371E-7,
      -2.1959984775055119E-7, -1.6153570148169649E-7, -1.1380671582905677E-7,
      -7.65909207167681E-8, -4.93053514682984E-8, -3.0883139823245642E-8,
      -1.9957835162287552E-8, -1.4964855798993426E-8, -1.4167336748194259E-8,
      -1.5598777158842836E-8, -1.6847016226506639E-8, -1.4293816796850933E-8,
      -1.4293816796850933E-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter_dag_3_DW.FromWorkspace_PWORK_e.TimePtr = (void *) pTimeValues0;
    helicopter_dag_3_DW.FromWorkspace_PWORK_e.DataPtr = (void *) pDataValues0;
    helicopter_dag_3_DW.FromWorkspace_IWORK_p.PrevIndex = 0;
  }

  /* Start for If: '<S3>/If' */
  helicopter_dag_3_DW.If_ActiveSubsystem = -1;

  /* Start for S-Function (game_controller_block): '<S6>/Game Controller' */

  /* S-Function Block: helicopter_dag_3/Joystick/Game Controller (game_controller_block) */
  {
    if (helicopter_dag_3_P.GameController_Enabled) {
      t_double deadzone[6];
      t_double saturation[6];
      t_int index;
      t_error result;
      for (index = 0; index < 6; index++) {
        deadzone[index] = -1;
      }

      for (index = 0; index < 6; index++) {
        saturation[index] = -1;
      }

      result = game_controller_open
        (helicopter_dag_3_P.GameController_ControllerNumber,
         helicopter_dag_3_P.GameController_BufferSize, deadzone, saturation,
         helicopter_dag_3_P.GameController_AutoCenter, 0, 1.0,
         &helicopter_dag_3_DW.GameController_Controller);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
      }
    }
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helicopter_dag_3_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helicopter_dag_3_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helicopter_dag_3_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  helicopter_dag_3_X.Integrator_CSTATE = helicopter_dag_3_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  helicopter_dag_3_DW.TimeStampA = (rtInf);
  helicopter_dag_3_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter_dag_3_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter_dag_3/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter_dag_3_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter_dag_3_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter_dag_3_P.HILInitialize_AOTerminate && !is_switching) ||
        (helicopter_dag_3_P.HILInitialize_AOExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_dag_3_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_dag_3_P.HILInitialize_AOFinal;
        }
      }

      num_final_analog_outputs = 8U;
    } else {
      num_final_analog_outputs = 0;
    }

    if ((helicopter_dag_3_P.HILInitialize_POTerminate && !is_switching) ||
        (helicopter_dag_3_P.HILInitialize_POExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_dag_3_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_dag_3_P.HILInitialize_POFinal;
        }
      }

      num_final_pwm_outputs = 8U;
    } else {
      num_final_pwm_outputs = 0;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter_dag_3_DW.HILInitialize_Card
                         , helicopter_dag_3_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , helicopter_dag_3_P.HILInitialize_POChannels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter_dag_3_DW.HILInitialize_AOVoltages[0]
                         , &helicopter_dag_3_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helicopter_dag_3_DW.HILInitialize_Card,
            helicopter_dag_3_P.HILInitialize_AOChannels,
            num_final_analog_outputs,
            &helicopter_dag_3_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter_dag_3_DW.HILInitialize_Card,
            helicopter_dag_3_P.HILInitialize_POChannels, num_final_pwm_outputs,
            &helicopter_dag_3_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_dag_3_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter_dag_3_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter_dag_3_DW.HILInitialize_Card);
    hil_close(helicopter_dag_3_DW.HILInitialize_Card);
    helicopter_dag_3_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for S-Function (game_controller_block): '<S6>/Game Controller' */

  /* S-Function Block: helicopter_dag_3/Joystick/Game Controller (game_controller_block) */
  {
    if (helicopter_dag_3_P.GameController_Enabled) {
      game_controller_close(helicopter_dag_3_DW.GameController_Controller);
      helicopter_dag_3_DW.GameController_Controller = NULL;
    }
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  helicopter_dag_3_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter_dag_3_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  helicopter_dag_3_initialize();
}

void MdlTerminate(void)
{
  helicopter_dag_3_terminate();
}

/* Registration function */
RT_MODEL_helicopter_dag_3_T *helicopter_dag_3(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter_dag_3_P.Integrator_UpperSat = rtInf;
  helicopter_dag_3_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter_dag_3_M, 0,
                sizeof(RT_MODEL_helicopter_dag_3_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter_dag_3_M->solverInfo,
                          &helicopter_dag_3_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter_dag_3_M->solverInfo, &rtmGetTPtr(helicopter_dag_3_M));
    rtsiSetStepSizePtr(&helicopter_dag_3_M->solverInfo,
                       &helicopter_dag_3_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter_dag_3_M->solverInfo, &helicopter_dag_3_M->derivs);
    rtsiSetContStatesPtr(&helicopter_dag_3_M->solverInfo, (real_T **)
                         &helicopter_dag_3_M->contStates);
    rtsiSetNumContStatesPtr(&helicopter_dag_3_M->solverInfo,
      &helicopter_dag_3_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&helicopter_dag_3_M->solverInfo,
      &helicopter_dag_3_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&helicopter_dag_3_M->solverInfo,
      &helicopter_dag_3_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&helicopter_dag_3_M->solverInfo,
      &helicopter_dag_3_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&helicopter_dag_3_M->solverInfo, (&rtmGetErrorStatus
      (helicopter_dag_3_M)));
    rtsiSetRTModelPtr(&helicopter_dag_3_M->solverInfo, helicopter_dag_3_M);
  }

  rtsiSetSimTimeStep(&helicopter_dag_3_M->solverInfo, MAJOR_TIME_STEP);
  helicopter_dag_3_M->intgData.y = helicopter_dag_3_M->odeY;
  helicopter_dag_3_M->intgData.f[0] = helicopter_dag_3_M->odeF[0];
  helicopter_dag_3_M->intgData.f[1] = helicopter_dag_3_M->odeF[1];
  helicopter_dag_3_M->intgData.f[2] = helicopter_dag_3_M->odeF[2];
  helicopter_dag_3_M->contStates = ((real_T *) &helicopter_dag_3_X);
  rtsiSetSolverData(&helicopter_dag_3_M->solverInfo, (void *)
                    &helicopter_dag_3_M->intgData);
  rtsiSetSolverName(&helicopter_dag_3_M->solverInfo,"ode3");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter_dag_3_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter_dag_3_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter_dag_3_M->Timing.sampleTimes =
      (&helicopter_dag_3_M->Timing.sampleTimesArray[0]);
    helicopter_dag_3_M->Timing.offsetTimes =
      (&helicopter_dag_3_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter_dag_3_M->Timing.sampleTimes[0] = (0.0);
    helicopter_dag_3_M->Timing.sampleTimes[1] = (0.02);

    /* task offsets */
    helicopter_dag_3_M->Timing.offsetTimes[0] = (0.0);
    helicopter_dag_3_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter_dag_3_M, &helicopter_dag_3_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter_dag_3_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter_dag_3_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter_dag_3_M, 25.0);
  helicopter_dag_3_M->Timing.stepSize0 = 0.02;
  helicopter_dag_3_M->Timing.stepSize1 = 0.02;

  /* External mode info */
  helicopter_dag_3_M->Sizes.checksums[0] = (515554214U);
  helicopter_dag_3_M->Sizes.checksums[1] = (3590486005U);
  helicopter_dag_3_M->Sizes.checksums[2] = (3118998705U);
  helicopter_dag_3_M->Sizes.checksums[3] = (3295601081U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[2];
    helicopter_dag_3_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = (sysRanDType *)
      &helicopter_dag_3_DW.IfActionSubsystem_SubsysRanBC;
    rteiSetModelMappingInfoPtr(helicopter_dag_3_M->extModeInfo,
      &helicopter_dag_3_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter_dag_3_M->extModeInfo,
                        helicopter_dag_3_M->Sizes.checksums);
    rteiSetTPtr(helicopter_dag_3_M->extModeInfo, rtmGetTPtr(helicopter_dag_3_M));
  }

  helicopter_dag_3_M->solverInfoPtr = (&helicopter_dag_3_M->solverInfo);
  helicopter_dag_3_M->Timing.stepSize = (0.02);
  rtsiSetFixedStepSize(&helicopter_dag_3_M->solverInfo, 0.02);
  rtsiSetSolverMode(&helicopter_dag_3_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter_dag_3_M->blockIO = ((void *) &helicopter_dag_3_B);

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_3_B.TmpSignalConversionAtToWorkspac[i] = 0.0;
    }

    helicopter_dag_3_B.FromWorkspace[0] = 0.0;
    helicopter_dag_3_B.FromWorkspace[1] = 0.0;
    helicopter_dag_3_B.FromWorkspace[2] = 0.0;
    helicopter_dag_3_B.FromWorkspace[3] = 0.0;
    helicopter_dag_3_B.TravelCounttorad = 0.0;
    helicopter_dag_3_B.Gain = 0.0;
    helicopter_dag_3_B.Gain_d = 0.0;
    helicopter_dag_3_B.PitchCounttorad = 0.0;
    helicopter_dag_3_B.Gain_i = 0.0;
    helicopter_dag_3_B.Gain_b = 0.0;
    helicopter_dag_3_B.ElevationCounttorad = 0.0;
    helicopter_dag_3_B.Gain_e = 0.0;
    helicopter_dag_3_B.Sum = 0.0;
    helicopter_dag_3_B.Gain_dg = 0.0;
    helicopter_dag_3_B.Subtract[0] = 0.0;
    helicopter_dag_3_B.Subtract[1] = 0.0;
    helicopter_dag_3_B.Subtract[2] = 0.0;
    helicopter_dag_3_B.Subtract[3] = 0.0;
    helicopter_dag_3_B.FromWorkspace1[0] = 0.0;
    helicopter_dag_3_B.FromWorkspace1[1] = 0.0;
    helicopter_dag_3_B.FromWorkspace1[2] = 0.0;
    helicopter_dag_3_B.FromWorkspace1[3] = 0.0;
    helicopter_dag_3_B.Subtract_a = 0.0;
    helicopter_dag_3_B.Gain_l = 0.0;
    helicopter_dag_3_B.BackmotorSaturation = 0.0;
    helicopter_dag_3_B.FrontmotorSaturation = 0.0;
    helicopter_dag_3_B.Joystick_gain_x = 0.0;
    helicopter_dag_3_B.Joystick_gain_y = 0.0;
    helicopter_dag_3_B.In1 = 0.0;
  }

  /* parameters */
  helicopter_dag_3_M->defaultParam = ((real_T *)&helicopter_dag_3_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter_dag_3_X;
    helicopter_dag_3_M->contStates = (x);
    (void) memset((void *)&helicopter_dag_3_X, 0,
                  sizeof(X_helicopter_dag_3_T));
  }

  /* states (dwork) */
  helicopter_dag_3_M->dwork = ((void *) &helicopter_dag_3_DW);
  (void) memset((void *)&helicopter_dag_3_DW, 0,
                sizeof(DW_helicopter_dag_3_T));

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_3_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_3_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_3_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_3_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_3_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_3_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_3_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_3_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter_dag_3_DW.TimeStampA = 0.0;
  helicopter_dag_3_DW.LastUAtTimeA = 0.0;
  helicopter_dag_3_DW.TimeStampB = 0.0;
  helicopter_dag_3_DW.LastUAtTimeB = 0.0;
  helicopter_dag_3_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter_dag_3_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter_dag_3_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 17;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter_dag_3_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter_dag_3_M->Sizes.numPeriodicContStates = (0);
                                      /* Number of periodic continuous states */
  helicopter_dag_3_M->Sizes.numY = (0);/* Number of model outputs */
  helicopter_dag_3_M->Sizes.numU = (0);/* Number of model inputs */
  helicopter_dag_3_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter_dag_3_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter_dag_3_M->Sizes.numBlocks = (79);/* Number of blocks */
  helicopter_dag_3_M->Sizes.numBlockIO = (21);/* Number of block outputs */
  helicopter_dag_3_M->Sizes.numBlockPrms = (164);/* Sum of parameter "widths" */
  return helicopter_dag_3_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
