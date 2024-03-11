/*
 * helicopter_dag_3.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "helicopter_dag_3".
 *
 * Model version              : 11.10
 * Simulink Coder version : 9.4 (R2020b) 29-Jul-2020
 * C source code generated on : Sun Mar 10 10:50:15 2024
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
 * Writes out MAT-file header.  Returns success or failure.
 * Returns:
 *      0 - success
 *      1 - failure
 */
int_T rt_WriteMat4FileHeader(FILE *fp, int32_T m, int32_T n, const char *name)
{
  typedef enum { ELITTLE_ENDIAN, EBIG_ENDIAN } ByteOrder;

  int16_T one = 1;
  ByteOrder byteOrder = (*((int8_T *)&one)==1) ? ELITTLE_ENDIAN : EBIG_ENDIAN;
  int32_T type = (byteOrder == ELITTLE_ENDIAN) ? 0: 1000;
  int32_T imagf = 0;
  int32_T name_len = (int32_T)strlen(name) + 1;
  if ((fwrite(&type, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&m, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&n, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&imagf, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&name_len, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(name, sizeof(char), name_len, fp) == 0)) {
    return(1);
  } else {
    return(0);
  }
}                                      /* end rt_WriteMat4FileHeader */

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
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_DeadZoney;
  real_T rtb_DeadZonex;
  real_T lastTime;
  real_T rtb_Clock;
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
    /* ToFile: '<Root>/To File' */
    {
      if (!(++helicopter_dag_3_DW.ToFile_IWORK.Decimation % 1) &&
          (helicopter_dag_3_DW.ToFile_IWORK.Count * (4 + 1)) + 1 < 100000000 ) {
        FILE *fp = (FILE *) helicopter_dag_3_DW.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[4 + 1];
          helicopter_dag_3_DW.ToFile_IWORK.Decimation = 0;
          u[0] = helicopter_dag_3_M->Timing.t[1];
          u[1] = helicopter_dag_3_B.Subtract[0];
          u[2] = helicopter_dag_3_B.Subtract[1];
          u[3] = helicopter_dag_3_B.Subtract[2];
          u[4] = helicopter_dag_3_B.Subtract[3];
          if (fwrite(u, sizeof(real_T), 4 + 1, fp) != 4 + 1) {
            rtmSetErrorStatus(helicopter_dag_3_M,
                              "Error writing to MAT-file plots/x_k_T5.mat");
            return;
          }

          if (((++helicopter_dag_3_DW.ToFile_IWORK.Count) * (4 + 1))+1 >=
              100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file plots/x_k_T5.mat.\n");
          }
        }
      }
    }
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
    } else if (t >= pTimeValues[116]) {
      currTimeIndex = 115;
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
              pDataValues += 117;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter_dag_3_B.FromWorkspace1[0])[elIdx] =
                pDataValues[currTimeIndex + 1];
              pDataValues += 117;
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
            pDataValues += 117;
          }
        }
      }
    }
  }

  /* FromWorkspace: '<S5>/From Workspace' */
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
    } else if (t >= pTimeValues[116]) {
      currTimeIndex = 115;
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
          helicopter_dag_3_B.FromWorkspace = pDataValues[currTimeIndex];
        } else {
          helicopter_dag_3_B.FromWorkspace = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        helicopter_dag_3_B.FromWorkspace = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 117;
      }
    }
  }

  /* Step: '<S5>/Step' */
  if (helicopter_dag_3_M->Timing.t[0] < helicopter_dag_3_P.T_padding) {
    rtb_Clock = helicopter_dag_3_P.Step_Y0;
  } else {
    rtb_Clock = helicopter_dag_3_P.Step_YFinal;
  }

  /* End of Step: '<S5>/Step' */

  /* Switch: '<S5>/Switch' incorporates:
   *  Constant: '<S5>/Constant'
   *  Gain: '<S5>/Gain'
   *  Sum: '<S5>/Subtract'
   *  Sum: '<S5>/Subtract1'
   */
  if (rtb_Clock > helicopter_dag_3_P.Switch_Threshold) {
    rtb_Clock = helicopter_dag_3_B.FromWorkspace -
      ((((helicopter_dag_3_B.Subtract[0] - helicopter_dag_3_B.FromWorkspace1[0])
         * helicopter_dag_3_P.K1[0] + (helicopter_dag_3_B.Subtract[1] -
          helicopter_dag_3_B.FromWorkspace1[1]) * helicopter_dag_3_P.K1[1]) +
        (helicopter_dag_3_B.Subtract[2] - helicopter_dag_3_B.FromWorkspace1[2]) *
        helicopter_dag_3_P.K1[2]) + (helicopter_dag_3_B.Subtract[3] -
        helicopter_dag_3_B.FromWorkspace1[3]) * helicopter_dag_3_P.K1[3]);
  } else {
    rtb_Clock = helicopter_dag_3_P.Constant_Value;
  }

  /* End of Switch: '<S5>/Switch' */

  /* Saturate: '<S5>/Saturation' */
  if (rtb_Clock > helicopter_dag_3_P.uu) {
    rtb_Clock = helicopter_dag_3_P.uu;
  } else {
    if (rtb_Clock < helicopter_dag_3_P.ul) {
      rtb_Clock = helicopter_dag_3_P.ul;
    }
  }

  /* End of Saturate: '<S5>/Saturation' */

  /* Gain: '<S5>/Gain1' */
  helicopter_dag_3_B.Gain1 = helicopter_dag_3_P.Gain1_Gain_g * rtb_Clock;
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
    /* ToFile: '<Root>/To File1' */
    {
      if (!(++helicopter_dag_3_DW.ToFile1_IWORK.Decimation % 1) &&
          (helicopter_dag_3_DW.ToFile1_IWORK.Count * (1 + 1)) + 1 < 100000000 )
      {
        FILE *fp = (FILE *) helicopter_dag_3_DW.ToFile1_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[1 + 1];
          helicopter_dag_3_DW.ToFile1_IWORK.Decimation = 0;
          u[0] = helicopter_dag_3_M->Timing.t[1];
          u[1] = helicopter_dag_3_B.Gain1;
          if (fwrite(u, sizeof(real_T), 1 + 1, fp) != 1 + 1) {
            rtmSetErrorStatus(helicopter_dag_3_M,
                              "Error writing to MAT-file plots/u_k_T5.mat");
            return;
          }

          if (((++helicopter_dag_3_DW.ToFile1_IWORK.Count) * (1 + 1))+1 >=
              100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file plots/u_k_T5.mat.\n");
          }
        }
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helicopter_dag_3_DW.FromWorkspace_PWORK_a.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter_dag_3_DW.FromWorkspace_PWORK_a.TimePtr;
    int_T currTimeIndex = helicopter_dag_3_DW.FromWorkspace_IWORK_d.PrevIndex;
    real_T t = helicopter_dag_3_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[116]) {
      currTimeIndex = 115;
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

    helicopter_dag_3_DW.FromWorkspace_IWORK_d.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter_dag_3_B.FromWorkspace_j[0])[elIdx] =
                pDataValues[currTimeIndex];
              pDataValues += 117;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter_dag_3_B.FromWorkspace_j[0])[elIdx] =
                pDataValues[currTimeIndex + 1];
              pDataValues += 117;
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
            (&helicopter_dag_3_B.FromWorkspace_j[0])[elIdx] = (real_T)
              rtInterpolate(d1, d2, f1, f2);
            pDataValues += 117;
          }
        }
      }
    }
  }

  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
  }

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S7>/K_pd'
   *  Gain: '<S7>/K_pp'
   *  Sum: '<S7>/Sum2'
   *  Sum: '<S7>/Sum3'
   */
  rtb_Clock = ((helicopter_dag_3_B.Gain1 - rtb_Gain1_idx_2) *
               helicopter_dag_3_P.K_pp - helicopter_dag_3_P.K_pd *
               rtb_Gain1_idx_3) + helicopter_dag_3_P.Vd_ff;

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
  rtb_Gain1_idx_2 = ((helicopter_dag_3_P.K_ep * rtb_Gain1_idx_3 +
                      helicopter_dag_3_X.Integrator_CSTATE) -
                     helicopter_dag_3_P.Gain1_Gain * helicopter_dag_3_B.Gain_dg *
                     helicopter_dag_3_P.K_ed) + helicopter_dag_3_P.Vs_ff;

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
    lastTime = helicopter_dag_3_DW.TimeStampA;
    lastU = &helicopter_dag_3_DW.LastUAtTimeA;
    if (helicopter_dag_3_DW.TimeStampA < helicopter_dag_3_DW.TimeStampB) {
      if (helicopter_dag_3_DW.TimeStampB < rtb_Gain1_idx_3) {
        lastTime = helicopter_dag_3_DW.TimeStampB;
        lastU = &helicopter_dag_3_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter_dag_3_DW.TimeStampA >= rtb_Gain1_idx_3) {
        lastTime = helicopter_dag_3_DW.TimeStampB;
        lastU = &helicopter_dag_3_DW.LastUAtTimeB;
      }
    }

    rtb_Gain1_idx_3 = (helicopter_dag_3_B.PitchCounttorad - *lastU) /
      (rtb_Gain1_idx_3 - lastTime);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S13>/Gain' */
  helicopter_dag_3_B.Gain_l = helicopter_dag_3_P.Gain_Gain_a1 * rtb_Gain1_idx_3;
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Gain1_idx_3 = (rtb_Gain1_idx_2 - rtb_Clock) *
    helicopter_dag_3_P.Backgain_Gain;

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Gain1_idx_3 > helicopter_dag_3_P.BackmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Back motor: Saturation' */
    helicopter_dag_3_B.BackmotorSaturation =
      helicopter_dag_3_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Gain1_idx_3 < helicopter_dag_3_P.BackmotorSaturation_LowerSat)
  {
    /* Saturate: '<S4>/Back motor: Saturation' */
    helicopter_dag_3_B.BackmotorSaturation =
      helicopter_dag_3_P.BackmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Back motor: Saturation' */
    helicopter_dag_3_B.BackmotorSaturation = rtb_Gain1_idx_3;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_dag_3_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Gain1_idx_3 = (rtb_Clock + rtb_Gain1_idx_2) *
    helicopter_dag_3_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Gain1_idx_3 > helicopter_dag_3_P.FrontmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Front motor: Saturation' */
    helicopter_dag_3_B.FrontmotorSaturation =
      helicopter_dag_3_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Gain1_idx_3 < helicopter_dag_3_P.FrontmotorSaturation_LowerSat)
  {
    /* Saturate: '<S4>/Front motor: Saturation' */
    helicopter_dag_3_B.FrontmotorSaturation =
      helicopter_dag_3_P.FrontmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Front motor: Saturation' */
    helicopter_dag_3_B.FrontmotorSaturation = rtb_Gain1_idx_3;
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

  /* Start for ToFile: '<Root>/To File' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "plots/x_k_T5.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helicopter_dag_3_M,
                        "Error creating .mat file plots/x_k_T5.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp, 4 + 1, 0, "ans")) {
      rtmSetErrorStatus(helicopter_dag_3_M,
                        "Error writing mat file header to file plots/x_k_T5.mat");
      return;
    }

    helicopter_dag_3_DW.ToFile_IWORK.Count = 0;
    helicopter_dag_3_DW.ToFile_IWORK.Decimation = -1;
    helicopter_dag_3_DW.ToFile_PWORK.FilePtr = fp;
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
      28.0, 28.25, 28.5, 28.75, 29.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1394479484589883, 3.1329867904146482,
      3.1206424903102814, 3.1015123173347732, 3.0752327569037745,
      3.0418250979055075, 3.0015635254379678, 2.9548768161827845,
      2.9022801718193505, 2.84433006112189, 2.7815952850372394,
      2.7146389546665808, 2.6440075882723075, 2.570224754851639,
      2.4937875776335625, 2.4151650167838437, 2.3347972500183078,
      2.2530957265356357, 2.1704436317073914, 2.0871966008687859,
      2.0036835827728603, 1.9202077913588038, 1.8370477076884313,
      1.754458108006232, 1.67267110245156, 1.5918971741843051,
      1.5123262119073435, 1.4341285307832294, 1.3574558780317016,
      1.2824424203476013, 1.2092057108682326, 1.1378476338476919,
      1.0684553255243381, 1.001102069932603, 0.93584816863334808,
      0.872741783530823, 0.8118197521167092, 0.75310837463743774,
      0.69662417282287525, 0.64237461994444223, 0.59035884209000766,
      0.54056829065239187, 0.49298738612867254, 0.44759413341927384,
      0.40436070889950176, 0.36325401961219189, 0.32423623499885523,
      0.28726529164851811, 0.2522953715987164, 0.21927735477216878,
      0.18815924617585744, 0.1588865785269182, 0.13140279100220195,
      0.10564958483592643, 0.0815672565127931, 0.05909500932258617,
      0.038171244056881506, 0.01873382963934225, 0.000720354488425691,
      -0.015931640584578844, -0.031284447138177837, -0.045399990721664181,
      -0.058339648556432346, -0.070164082486200888, -0.080933086413504274,
      -0.090705447450565357, -0.099538820025200828, -0.10748961219657842,
      -0.11461288345127646, -0.12096225326704785, -0.1265898197498227,
      -0.13154608766867162, -0.13587990523357973, -0.13963840898184846,
      -0.14286697616066896, -0.14560918401583378, -0.14790677541965216,
      -0.1497996302949352, -0.15132574231653298, -0.15252120039758421,
      -0.15342017449484358, -0.15405490529701452, -0.15445569739335441,
      -0.15465091555932572, -0.154666983845715, -0.15452838722394591,
      -0.15425767563385287, -0.15387547041790162, -0.15340047333444723,
      -0.15284947866332779, -0.15223738941165424, -0.15157723938460996,
      -0.15088022402347906, -0.1501557445677689, -0.14941147237703908,
      -0.1486534430856829, -0.14788619306521039, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -0.008578820523219283, -0.02584463217736142, -0.049377200417466877,
      -0.0765206919020325, -0.1051182417239961, -0.13363063599306871,
      -0.1610462898701584, -0.18674683702073208, -0.21038657745373604,
      -0.23180044278984285, -0.25093910433860084, -0.2678253214826356,
      -0.28252546557709268, -0.29513133368267391, -0.30574870887230421,
      -0.31449024339887621, -0.32147106706214179, -0.32680609393068738,
      -0.3306083793129771, -0.33298812335442246, -0.33405207238370305,
      -0.33390316565622508, -0.33264033468149018, -0.33035839872879641,
      -0.32714802221868827, -0.32309571306901952, -0.31828384910784607,
      -0.31279072449645545, -0.30669061100611178, -0.30005383073640091,
      -0.29294683791747445, -0.28543230808216236, -0.27756923329341493,
      -0.26941302236694048, -0.26101560519702005, -0.25242554041010029,
      -0.2436881256564549, -0.23484550991708614, -0.22593680725825005,
      -0.21699821151373189, -0.20806311141773826, -0.19916220575046312,
      -0.19032361809487719, -0.18157301083759481, -0.17293369807908848,
      -0.16442675714923954, -0.15607113845334666, -0.14788377340134853,
      -0.13987968019920674, -0.13207206730619045, -0.12447243438524541,
      -0.11709067059575691, -0.10993515009886497, -0.10301282466510213,
      -0.096329313292533345, -0.0898889887608277, -0.083695061062818682,
      -0.077749657670157027, -0.072053900603666232, -0.066607980292018129,
      -0.061411226214395973, -0.0564621743339454, -0.051758631339072665,
      -0.047297735719074134, -0.043076015709213579, -0.039089444148244341,
      -0.035333490298541888, -0.031803168685510343, -0.028493085018792188,
      -0.02539747926308554, -0.022510265931099314, -0.0198250716753958,
      -0.017335270259632503, -0.015034014993074897, -0.012914268715281924,
      -0.010968831420659273, -0.0091903656152736261, -0.0075714195011321048,
      -0.0061044480863911245, -0.0047818323242049725, -0.0035958963890374059,
      -0.002538923208683784, -0.0016031683853595527, -0.00078087266388525412,
      -6.4273145557180756E-5, 0.000554386487076434, 0.0010828463603721648,
      0.0015288208638049663, 0.0018999883338175678, 0.002203978684477758,
      0.0024483570066941907, 0.0026406001081770059, 0.0027880614445235973,
      0.0028979178228406371, 0.0029770887629192323, 0.003032117165424804,
      0.0030690000818899307, 0.0030929663720603027, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.060631827007488148, 0.12202816255732984, 0.16631920468713923,
      0.19183983108359426, 0.20211656007369805, 0.20151471314874658,
      0.19376337092882057, 0.18164165162486834, 0.16707665681499839,
      0.15134502174348685, 0.13526475032763452, 0.11934533353584398,
      0.10389500413346109, 0.089093461296690935, 0.075039552818343291,
      0.061781827438931236, 0.049337795513422433, 0.037705932336332904,
      0.026873100901022351, 0.016819122005487075, 0.0075195853921584677,
      -0.0010524158788496329, -0.0089252070249395965, -0.016127851789310421,
      -0.022689715056118831, -0.028640173398918611, -0.034008416715023526,
      -0.038823306801505919, -0.043113272374805611, -0.046906228205875033,
      -0.0502295109186458, -0.053109826902614232, -0.055573209516088629,
      -0.057644983782049242, -0.059349737388463053, -0.0607112971688184,
      -0.061752710457263105, -0.062496230845657974, -0.062963307952324521,
      -0.063174580865062446, -0.063149874956360319, -0.062908201793824,
      -0.062467761887899487, -0.061845950034446318, -0.061059363023048929,
      -0.060123809493971947, -0.059054321737843529, -0.057865169242766967,
      -0.056569873803755166, -0.055181226019246132, -0.05371130300902438,
      -0.052171487197167354, -0.050572486012658957, -0.048924352369069113,
      -0.047236505793181782, -0.045517754080659611, -0.043776315364761065,
      -0.042019840491765825, -0.040255435604115986, -0.03848968483934273,
      -0.0367286730596148, -0.0349780085332238, -0.033242845495500772,
      -0.031527906522553772, -0.029837504656814784, -0.028175565228701016,
      -0.026545647323723498, -0.024950964849121071, -0.023394407158561004,
      -0.021878559197623493, -0.020405721136670607, -0.018977927461278021,
      -0.017596965493640981, -0.016264393321205217, -0.014981557111112043,
      -0.013749607790710217, -0.012569517075069925, -0.011442092821628558,
      -0.010367993688930377, -0.0093477430694466346, -0.0083817422532773689,
      -0.0074702827561213825, -0.0066135577047414795, -0.0058116721055302412,
      -0.005064651709525636, -0.00437245000295261, -0.0037349525533338923,
      -0.003151977462219735, -0.0026232699205479548, -0.0021484876974660061,
      -0.0017271726476426608, -0.0013587008191482397, -0.001042200406360716,
      -0.00077642292521418188, -0.00055954996722584749, -0.00038891973225507748,
      -0.00026067436710353054, -0.0001693845856115539, -0.00010785488567366386,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.24252730802984834, 0.245585342199382, 0.17716416851924846,
      0.10208250558582616, 0.041106915960418017, -0.0024073876998062316,
      -0.031005368879705734, -0.048486877215811894, -0.058259979239483471,
      -0.062926540286050034, -0.064321085663413047, -0.063677667167166288,
      -0.06180131760953541, -0.059206171347084259, -0.056215633913394084,
      -0.053030901517651517, -0.049776127702038125, -0.0465274527083607,
      -0.043331325741245111, -0.040215915582143394, -0.037198146453317017,
      -0.03428800508403413, -0.03149116458436186, -0.028810579057484836,
      -0.026247453067235485, -0.023801833371200506, -0.021472973264421182,
      -0.019259560345930564, -0.017159862293199952, -0.015171823324278806,
      -0.013293130851084013, -0.011521263935874479, -0.0098535304538980989,
      -0.0082870970638430838, -0.0068190144256558855, -0.0054462391214215589,
      -0.0041656531537790937, -0.002974081553579571, -0.0018683084266665339,
      -0.00084509165095172825, 9.8823634808343789E-5, 0.00096669265014522553,
      0.0017617596236980714, 0.0024872474138129624, 0.003146348045589584,
      0.0037422141163081591, 0.0042779510245139363, 0.004756609980306554,
      0.0051811817560474633, 0.0055545911380365026, 0.0058796920408872525,
      0.0061592632474284925, 0.0063960047380341445, 0.0065925345743597661,
      0.0067513863035496592, 0.0068750068500890828, 0.0069657548635946856,
      0.0070258994919814889, 0.0070576195505996225, 0.0070630030590934541,
      0.0070440471189120021, 0.0070026581055644584, 0.0069406521508924065,
      0.0068597558917883978, 0.006761607462956302, 0.0066477577124553584,
      0.0065196716199105043, 0.0063787298984101333, 0.0062262307622403396,
      0.0060633918437504536, 0.0058913522438119315, 0.0057111747015706274,
      0.005523847870548413, 0.0053302886897431912, 0.005131344840372685,
      0.0049277972816077677, 0.0047203628625616676, 0.0045096970137651895,
      0.0042963965307933616, 0.0040810024779348855, 0.00386400326467753,
      0.0036458379886239551, 0.0034269002055196032, 0.0032075423968451217,
      0.0029880815840184886, 0.0027688068262927109, 0.0025499897984746734,
      0.0023319003644566546, 0.0021148301666872576, 0.001899128892327665,
      0.0016852601992936604, 0.0014738873139776179, 0.0012660016511499925,
      0.0010631099245865328, 0.0008674918319534642, 0.00068252093988297726,
      0.00051298146060631635, 0.00036515912596824033, 0.00024611879975124388,
      0.00016075865679742276, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

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
      28.0, 28.25, 28.5, 28.75, 29.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.29941642966649074, 0.27325013752642935, 0.24903384466094602,
      0.22618382766090528, 0.20446570518323004, 0.18379281960723914,
      0.16413553524687929, 0.14548219422548903, 0.12782350485734029,
      0.11114715488153726, 0.09543661775968415, 0.08067149260667017,
      0.066828267743911818, 0.053881087133217265, 0.041802390826445657,
      0.030563414227173347, 0.020134568188379887, 0.010485728749851764,
      0.001586461599340061, -0.006593799772329656, -0.014085608846015618,
      -0.020919384708957933, -0.027125308528775816, -0.032733233902364489,
      -0.03777260845887237, -0.0422724049553298, -0.046261060645378105,
      -0.049766424040261104, -0.052815708392314864, -0.0554354513631824,
      -0.057651480422528212, -0.0594888835771904, -0.06097198506717677,
      -0.0621243256908276, -0.062968647441114745, -0.063526882151006769,
      -0.063820143859619427, -0.0638687246233881, -0.063692093508261216,
      -0.063308898510207379, -0.062736971162320188, -0.061993333597575639,
      -0.061094207846874471, -0.060055027162412788, -0.058890449166654923,
      -0.0576143706372253, -0.05623994374788488, -0.054779593595385645,
      -0.05324503685140336, -0.051647301387913291, -0.049996746733281916,
      -0.048303085224999875, -0.046575403733356446, -0.044822185838456319,
      -0.043051334350790615, -0.041270194073101396, -0.039485574708508442,
      -0.037703773826805675, -0.03593059980747576, -0.034171394684317513,
      -0.0324310568226317, -0.030714063365672395, -0.029024492392533796,
      -0.027366044734830286, -0.025742065404417924, -0.024155564589023626,
      -0.022609238176981794, -0.021105487776332144, -0.019646440197304882,
      -0.018233966370693744, -0.016869699677775074, -0.015555053670229468,
      -0.014291239160879154, -0.013079280667835924, -0.011920032195618546,
      -0.010814192336550454, -0.00976231867358035, -0.0087648414604524216,
      -0.0078220765449259488, -0.0069342374823010422, -0.0061014467545058571,
      -0.0053237459556738642, -0.0046011047141945038, -0.003933427970574388,
      -0.0033205609833416583, -0.0027622910340601159, -0.0022583441590993125,
      -0.0018083742194845431, -0.0014119400468629584, -0.0010684640354804431,
      -0.00077716213034939585, -0.00053693054432002452, -0.00034616904301032481,
      -0.00020251593091191733, -0.00010247097488580437, -4.0902237945505959E-5,
      -1.0501366462967177E-5, -1.43704101129849E-6, -1.9053706168392637E-6,
      -1.3021450903094589E-6, -1.3021450903094589E-6, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0 } ;

    helicopter_dag_3_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter_dag_3_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helicopter_dag_3_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File1' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "plots/u_k_T5.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helicopter_dag_3_M,
                        "Error creating .mat file plots/u_k_T5.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp, 1 + 1, 0, "ans")) {
      rtmSetErrorStatus(helicopter_dag_3_M,
                        "Error writing mat file header to file plots/u_k_T5.mat");
      return;
    }

    helicopter_dag_3_DW.ToFile1_IWORK.Count = 0;
    helicopter_dag_3_DW.ToFile1_IWORK.Decimation = -1;
    helicopter_dag_3_DW.ToFile1_PWORK.FilePtr = fp;
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
      28.0, 28.25, 28.5, 28.75, 29.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1394479484589883, 3.1329867904146482,
      3.1206424903102814, 3.1015123173347732, 3.0752327569037745,
      3.0418250979055075, 3.0015635254379678, 2.9548768161827845,
      2.9022801718193505, 2.84433006112189, 2.7815952850372394,
      2.7146389546665808, 2.6440075882723075, 2.570224754851639,
      2.4937875776335625, 2.4151650167838437, 2.3347972500183078,
      2.2530957265356357, 2.1704436317073914, 2.0871966008687859,
      2.0036835827728603, 1.9202077913588038, 1.8370477076884313,
      1.754458108006232, 1.67267110245156, 1.5918971741843051,
      1.5123262119073435, 1.4341285307832294, 1.3574558780317016,
      1.2824424203476013, 1.2092057108682326, 1.1378476338476919,
      1.0684553255243381, 1.001102069932603, 0.93584816863334808,
      0.872741783530823, 0.8118197521167092, 0.75310837463743774,
      0.69662417282287525, 0.64237461994444223, 0.59035884209000766,
      0.54056829065239187, 0.49298738612867254, 0.44759413341927384,
      0.40436070889950176, 0.36325401961219189, 0.32423623499885523,
      0.28726529164851811, 0.2522953715987164, 0.21927735477216878,
      0.18815924617585744, 0.1588865785269182, 0.13140279100220195,
      0.10564958483592643, 0.0815672565127931, 0.05909500932258617,
      0.038171244056881506, 0.01873382963934225, 0.000720354488425691,
      -0.015931640584578844, -0.031284447138177837, -0.045399990721664181,
      -0.058339648556432346, -0.070164082486200888, -0.080933086413504274,
      -0.090705447450565357, -0.099538820025200828, -0.10748961219657842,
      -0.11461288345127646, -0.12096225326704785, -0.1265898197498227,
      -0.13154608766867162, -0.13587990523357973, -0.13963840898184846,
      -0.14286697616066896, -0.14560918401583378, -0.14790677541965216,
      -0.1497996302949352, -0.15132574231653298, -0.15252120039758421,
      -0.15342017449484358, -0.15405490529701452, -0.15445569739335441,
      -0.15465091555932572, -0.154666983845715, -0.15452838722394591,
      -0.15425767563385287, -0.15387547041790162, -0.15340047333444723,
      -0.15284947866332779, -0.15223738941165424, -0.15157723938460996,
      -0.15088022402347906, -0.1501557445677689, -0.14941147237703908,
      -0.1486534430856829, -0.14788619306521039, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -0.008578820523219283, -0.02584463217736142, -0.049377200417466877,
      -0.0765206919020325, -0.1051182417239961, -0.13363063599306871,
      -0.1610462898701584, -0.18674683702073208, -0.21038657745373604,
      -0.23180044278984285, -0.25093910433860084, -0.2678253214826356,
      -0.28252546557709268, -0.29513133368267391, -0.30574870887230421,
      -0.31449024339887621, -0.32147106706214179, -0.32680609393068738,
      -0.3306083793129771, -0.33298812335442246, -0.33405207238370305,
      -0.33390316565622508, -0.33264033468149018, -0.33035839872879641,
      -0.32714802221868827, -0.32309571306901952, -0.31828384910784607,
      -0.31279072449645545, -0.30669061100611178, -0.30005383073640091,
      -0.29294683791747445, -0.28543230808216236, -0.27756923329341493,
      -0.26941302236694048, -0.26101560519702005, -0.25242554041010029,
      -0.2436881256564549, -0.23484550991708614, -0.22593680725825005,
      -0.21699821151373189, -0.20806311141773826, -0.19916220575046312,
      -0.19032361809487719, -0.18157301083759481, -0.17293369807908848,
      -0.16442675714923954, -0.15607113845334666, -0.14788377340134853,
      -0.13987968019920674, -0.13207206730619045, -0.12447243438524541,
      -0.11709067059575691, -0.10993515009886497, -0.10301282466510213,
      -0.096329313292533345, -0.0898889887608277, -0.083695061062818682,
      -0.077749657670157027, -0.072053900603666232, -0.066607980292018129,
      -0.061411226214395973, -0.0564621743339454, -0.051758631339072665,
      -0.047297735719074134, -0.043076015709213579, -0.039089444148244341,
      -0.035333490298541888, -0.031803168685510343, -0.028493085018792188,
      -0.02539747926308554, -0.022510265931099314, -0.0198250716753958,
      -0.017335270259632503, -0.015034014993074897, -0.012914268715281924,
      -0.010968831420659273, -0.0091903656152736261, -0.0075714195011321048,
      -0.0061044480863911245, -0.0047818323242049725, -0.0035958963890374059,
      -0.002538923208683784, -0.0016031683853595527, -0.00078087266388525412,
      -6.4273145557180756E-5, 0.000554386487076434, 0.0010828463603721648,
      0.0015288208638049663, 0.0018999883338175678, 0.002203978684477758,
      0.0024483570066941907, 0.0026406001081770059, 0.0027880614445235973,
      0.0028979178228406371, 0.0029770887629192323, 0.003032117165424804,
      0.0030690000818899307, 0.0030929663720603027, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.060631827007488148, 0.12202816255732984, 0.16631920468713923,
      0.19183983108359426, 0.20211656007369805, 0.20151471314874658,
      0.19376337092882057, 0.18164165162486834, 0.16707665681499839,
      0.15134502174348685, 0.13526475032763452, 0.11934533353584398,
      0.10389500413346109, 0.089093461296690935, 0.075039552818343291,
      0.061781827438931236, 0.049337795513422433, 0.037705932336332904,
      0.026873100901022351, 0.016819122005487075, 0.0075195853921584677,
      -0.0010524158788496329, -0.0089252070249395965, -0.016127851789310421,
      -0.022689715056118831, -0.028640173398918611, -0.034008416715023526,
      -0.038823306801505919, -0.043113272374805611, -0.046906228205875033,
      -0.0502295109186458, -0.053109826902614232, -0.055573209516088629,
      -0.057644983782049242, -0.059349737388463053, -0.0607112971688184,
      -0.061752710457263105, -0.062496230845657974, -0.062963307952324521,
      -0.063174580865062446, -0.063149874956360319, -0.062908201793824,
      -0.062467761887899487, -0.061845950034446318, -0.061059363023048929,
      -0.060123809493971947, -0.059054321737843529, -0.057865169242766967,
      -0.056569873803755166, -0.055181226019246132, -0.05371130300902438,
      -0.052171487197167354, -0.050572486012658957, -0.048924352369069113,
      -0.047236505793181782, -0.045517754080659611, -0.043776315364761065,
      -0.042019840491765825, -0.040255435604115986, -0.03848968483934273,
      -0.0367286730596148, -0.0349780085332238, -0.033242845495500772,
      -0.031527906522553772, -0.029837504656814784, -0.028175565228701016,
      -0.026545647323723498, -0.024950964849121071, -0.023394407158561004,
      -0.021878559197623493, -0.020405721136670607, -0.018977927461278021,
      -0.017596965493640981, -0.016264393321205217, -0.014981557111112043,
      -0.013749607790710217, -0.012569517075069925, -0.011442092821628558,
      -0.010367993688930377, -0.0093477430694466346, -0.0083817422532773689,
      -0.0074702827561213825, -0.0066135577047414795, -0.0058116721055302412,
      -0.005064651709525636, -0.00437245000295261, -0.0037349525533338923,
      -0.003151977462219735, -0.0026232699205479548, -0.0021484876974660061,
      -0.0017271726476426608, -0.0013587008191482397, -0.001042200406360716,
      -0.00077642292521418188, -0.00055954996722584749, -0.00038891973225507748,
      -0.00026067436710353054, -0.0001693845856115539, -0.00010785488567366386,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.24252730802984834, 0.245585342199382, 0.17716416851924846,
      0.10208250558582616, 0.041106915960418017, -0.0024073876998062316,
      -0.031005368879705734, -0.048486877215811894, -0.058259979239483471,
      -0.062926540286050034, -0.064321085663413047, -0.063677667167166288,
      -0.06180131760953541, -0.059206171347084259, -0.056215633913394084,
      -0.053030901517651517, -0.049776127702038125, -0.0465274527083607,
      -0.043331325741245111, -0.040215915582143394, -0.037198146453317017,
      -0.03428800508403413, -0.03149116458436186, -0.028810579057484836,
      -0.026247453067235485, -0.023801833371200506, -0.021472973264421182,
      -0.019259560345930564, -0.017159862293199952, -0.015171823324278806,
      -0.013293130851084013, -0.011521263935874479, -0.0098535304538980989,
      -0.0082870970638430838, -0.0068190144256558855, -0.0054462391214215589,
      -0.0041656531537790937, -0.002974081553579571, -0.0018683084266665339,
      -0.00084509165095172825, 9.8823634808343789E-5, 0.00096669265014522553,
      0.0017617596236980714, 0.0024872474138129624, 0.003146348045589584,
      0.0037422141163081591, 0.0042779510245139363, 0.004756609980306554,
      0.0051811817560474633, 0.0055545911380365026, 0.0058796920408872525,
      0.0061592632474284925, 0.0063960047380341445, 0.0065925345743597661,
      0.0067513863035496592, 0.0068750068500890828, 0.0069657548635946856,
      0.0070258994919814889, 0.0070576195505996225, 0.0070630030590934541,
      0.0070440471189120021, 0.0070026581055644584, 0.0069406521508924065,
      0.0068597558917883978, 0.006761607462956302, 0.0066477577124553584,
      0.0065196716199105043, 0.0063787298984101333, 0.0062262307622403396,
      0.0060633918437504536, 0.0058913522438119315, 0.0057111747015706274,
      0.005523847870548413, 0.0053302886897431912, 0.005131344840372685,
      0.0049277972816077677, 0.0047203628625616676, 0.0045096970137651895,
      0.0042963965307933616, 0.0040810024779348855, 0.00386400326467753,
      0.0036458379886239551, 0.0034269002055196032, 0.0032075423968451217,
      0.0029880815840184886, 0.0027688068262927109, 0.0025499897984746734,
      0.0023319003644566546, 0.0021148301666872576, 0.001899128892327665,
      0.0016852601992936604, 0.0014738873139776179, 0.0012660016511499925,
      0.0010631099245865328, 0.0008674918319534642, 0.00068252093988297726,
      0.00051298146060631635, 0.00036515912596824033, 0.00024611879975124388,
      0.00016075865679742276, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter_dag_3_DW.FromWorkspace_PWORK_a.TimePtr = (void *) pTimeValues0;
    helicopter_dag_3_DW.FromWorkspace_PWORK_a.DataPtr = (void *) pDataValues0;
    helicopter_dag_3_DW.FromWorkspace_IWORK_d.PrevIndex = 0;
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

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helicopter_dag_3_DW.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "plots/x_k_T5.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_dag_3_M,
                          "Error closing MAT-file plots/x_k_T5.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helicopter_dag_3_M,
                          "Error reopening MAT-file plots/x_k_T5.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 4 + 1,
           helicopter_dag_3_DW.ToFile_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helicopter_dag_3_M,
                          "Error writing header for ans to MAT-file plots/x_k_T5.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_dag_3_M,
                          "Error closing MAT-file plots/x_k_T5.mat");
        return;
      }

      helicopter_dag_3_DW.ToFile_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for ToFile: '<Root>/To File1' */
  {
    FILE *fp = (FILE *) helicopter_dag_3_DW.ToFile1_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "plots/u_k_T5.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_dag_3_M,
                          "Error closing MAT-file plots/u_k_T5.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helicopter_dag_3_M,
                          "Error reopening MAT-file plots/u_k_T5.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 1 + 1,
           helicopter_dag_3_DW.ToFile1_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helicopter_dag_3_M,
                          "Error writing header for ans to MAT-file plots/u_k_T5.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_dag_3_M,
                          "Error closing MAT-file plots/u_k_T5.mat");
        return;
      }

      helicopter_dag_3_DW.ToFile1_PWORK.FilePtr = (NULL);
    }
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

  rtmSetTFinal(helicopter_dag_3_M, 20.0);
  helicopter_dag_3_M->Timing.stepSize0 = 0.02;
  helicopter_dag_3_M->Timing.stepSize1 = 0.02;

  /* External mode info */
  helicopter_dag_3_M->Sizes.checksums[0] = (2436993970U);
  helicopter_dag_3_M->Sizes.checksums[1] = (177235433U);
  helicopter_dag_3_M->Sizes.checksums[2] = (1580888293U);
  helicopter_dag_3_M->Sizes.checksums[3] = (1659806876U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[3];
    helicopter_dag_3_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = (sysRanDType *)
      &helicopter_dag_3_DW.IfActionSubsystem_SubsysRanBC;
    systemRan[2] = &rtAlwaysEnabled;
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
    helicopter_dag_3_B.FromWorkspace = 0.0;
    helicopter_dag_3_B.Gain1 = 0.0;
    helicopter_dag_3_B.FromWorkspace_j[0] = 0.0;
    helicopter_dag_3_B.FromWorkspace_j[1] = 0.0;
    helicopter_dag_3_B.FromWorkspace_j[2] = 0.0;
    helicopter_dag_3_B.FromWorkspace_j[3] = 0.0;
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
  helicopter_dag_3_M->Sizes.numBlocks = (85);/* Number of blocks */
  helicopter_dag_3_M->Sizes.numBlockIO = (21);/* Number of block outputs */
  helicopter_dag_3_M->Sizes.numBlockPrms = (172);/* Sum of parameter "widths" */
  return helicopter_dag_3_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
