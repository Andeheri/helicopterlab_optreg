/*
 * helicopter_dag_2.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "helicopter_dag_2".
 *
 * Model version              : 11.7
 * Simulink Coder version : 9.4 (R2020b) 29-Jul-2020
 * C source code generated on : Wed Feb 14 15:39:10 2024
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helicopter_dag_2.h"
#include "helicopter_dag_2_private.h"
#include "helicopter_dag_2_dt.h"

/* Block signals (default storage) */
B_helicopter_dag_2_T helicopter_dag_2_B;

/* Continuous states */
X_helicopter_dag_2_T helicopter_dag_2_X;

/* Block states (default storage) */
DW_helicopter_dag_2_T helicopter_dag_2_DW;

/* Real-time model */
static RT_MODEL_helicopter_dag_2_T helicopter_dag_2_M_;
RT_MODEL_helicopter_dag_2_T *const helicopter_dag_2_M = &helicopter_dag_2_M_;
static void rate_monotonic_scheduler(void);

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

time_T rt_SimUpdateDiscreteEvents(
  int_T rtmNumSampTimes, void *rtmTimingData, int_T *rtmSampleHitPtr, int_T
  *rtmPerTaskSampleHits )
{
  rtmSampleHitPtr[1] = rtmStepTask(helicopter_dag_2_M, 1);
  rtmSampleHitPtr[2] = rtmStepTask(helicopter_dag_2_M, 2);
  UNUSED_PARAMETER(rtmNumSampTimes);
  UNUSED_PARAMETER(rtmTimingData);
  UNUSED_PARAMETER(rtmPerTaskSampleHits);
  return(-1);
}

/*
 *   This function updates active task flag for each subrate
 * and rate transition flags for tasks that exchange data.
 * The function assumes rate-monotonic multitasking scheduler.
 * The function must be called at model base rate so that
 * the generated code self-manages all its subrates and rate
 * transition flags.
 */
static void rate_monotonic_scheduler(void)
{
  /* To ensure a deterministic data transfer between two rates,
   * data is transferred at the priority of a fast task and the frequency
   * of the slow task.  The following flags indicate when the data transfer
   * happens.  That is, a rate interaction flag is set true when both rates
   * will run, and false otherwise.
   */

  /* tid 1 shares data with slower tid rate: 2 */
  if (helicopter_dag_2_M->Timing.TaskCounters.TID[1] == 0) {
    helicopter_dag_2_M->Timing.RateInteraction.TID1_2 =
      (helicopter_dag_2_M->Timing.TaskCounters.TID[2] == 0);

    /* update PerTaskSampleHits matrix for non-inline sfcn */
    helicopter_dag_2_M->Timing.perTaskSampleHits[5] =
      helicopter_dag_2_M->Timing.RateInteraction.TID1_2;
  }

  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (helicopter_dag_2_M->Timing.TaskCounters.TID[2])++;
  if ((helicopter_dag_2_M->Timing.TaskCounters.TID[2]) > 4) {/* Sample time: [0.01s, 0.0s] */
    helicopter_dag_2_M->Timing.TaskCounters.TID[2] = 0;
  }
}

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  helicopter_dag_2_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function for TID0 */
void helicopter_dag_2_output0(void)    /* Sample time: [0.0s, 0.0s] */
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_DeadZonex;
  real_T lastTime;
  real_T rtb_Backgain;
  real_T rtb_Clock;
  real_T rtb_Derivative;
  real_T *lastU;
  int8_T rtAction;
  if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
    /* set solver stop time */
    if (!(helicopter_dag_2_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter_dag_2_M->solverInfo,
                            ((helicopter_dag_2_M->Timing.clockTickH0 + 1) *
        helicopter_dag_2_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter_dag_2_M->solverInfo,
                            ((helicopter_dag_2_M->Timing.clockTick0 + 1) *
        helicopter_dag_2_M->Timing.stepSize0 +
        helicopter_dag_2_M->Timing.clockTickH0 *
        helicopter_dag_2_M->Timing.stepSize0 * 4294967296.0));
    }

    {                                  /* Sample time: [0.0s, 0.0s] */
      rate_monotonic_scheduler();
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter_dag_2_M)) {
    helicopter_dag_2_M->Timing.t[0] = rtsiGetT(&helicopter_dag_2_M->solverInfo);
  }

  /* Reset subsysRan breadcrumbs */
  srClearBC(helicopter_dag_2_DW.IfActionSubsystem_SubsysRanBC);
  if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter_dag_2/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (helicopter_dag_2_DW.HILReadEncoderTimebase_Task, 1,
         &helicopter_dag_2_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter_dag_2_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter_dag_2_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_DeadZonex = helicopter_dag_2_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helicopter_dag_2_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter_dag_2_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helicopter_dag_2_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helicopter_dag_2_M->Timing.t[0];

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

    helicopter_dag_2_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          helicopter_dag_2_B.FromWorkspace = pDataValues[currTimeIndex];
        } else {
          helicopter_dag_2_B.FromWorkspace = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        helicopter_dag_2_B.FromWorkspace = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
    /* Gain: '<S4>/Travel: Count to rad' incorporates:
     *  Gain: '<S4>/Travel_gain'
     */
    helicopter_dag_2_B.TravelCounttorad = helicopter_dag_2_P.travel_gain *
      rtb_HILReadEncoderTimebase_o1 * helicopter_dag_2_P.TravelCounttorad_Gain;

    /* Gain: '<S13>/Gain' */
    helicopter_dag_2_B.Gain = helicopter_dag_2_P.Gain_Gain *
      helicopter_dag_2_B.TravelCounttorad;

    /* Gain: '<S4>/Pitch: Count to rad' */
    helicopter_dag_2_B.PitchCounttorad = helicopter_dag_2_P.PitchCounttorad_Gain
      * rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S10>/Gain' */
    helicopter_dag_2_B.Gain_i = helicopter_dag_2_P.Gain_Gain_a *
      helicopter_dag_2_B.PitchCounttorad;
  }

  /* Gain: '<S14>/Gain' incorporates:
   *  TransferFcn: '<S4>/Travel: Transfer Fcn'
   */
  helicopter_dag_2_B.Gain_d = (helicopter_dag_2_P.TravelTransferFcn_C *
    helicopter_dag_2_X.TravelTransferFcn_CSTATE +
    helicopter_dag_2_P.TravelTransferFcn_D * helicopter_dag_2_B.TravelCounttorad)
    * helicopter_dag_2_P.Gain_Gain_l;

  /* Gain: '<S11>/Gain' incorporates:
   *  TransferFcn: '<S4>/Pitch: Transfer Fcn'
   */
  helicopter_dag_2_B.Gain_b = (helicopter_dag_2_P.PitchTransferFcn_C *
    helicopter_dag_2_X.PitchTransferFcn_CSTATE +
    helicopter_dag_2_P.PitchTransferFcn_D * helicopter_dag_2_B.PitchCounttorad) *
    helicopter_dag_2_P.Gain_Gain_ae;
  if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' incorporates:
     *  Gain: '<S4>/Elevation_gain'
     */
    helicopter_dag_2_B.ElevationCounttorad = helicopter_dag_2_P.elevation_gain *
      rtb_DeadZonex * helicopter_dag_2_P.ElevationCounttorad_Gain;

    /* Gain: '<S8>/Gain' */
    helicopter_dag_2_B.Gain_e = helicopter_dag_2_P.Gain_Gain_lv *
      helicopter_dag_2_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter_dag_2_B.Sum = helicopter_dag_2_B.Gain_e +
      helicopter_dag_2_P.elavation_offsetdeg_Value;
  }

  /* Gain: '<S9>/Gain' incorporates:
   *  TransferFcn: '<S4>/Elevation: Transfer Fcn'
   */
  helicopter_dag_2_B.Gain_dg = (helicopter_dag_2_P.ElevationTransferFcn_C *
    helicopter_dag_2_X.ElevationTransferFcn_CSTATE +
    helicopter_dag_2_P.ElevationTransferFcn_D *
    helicopter_dag_2_B.ElevationCounttorad) * helicopter_dag_2_P.Gain_Gain_n;

  /* Gain: '<S2>/Gain1' */
  helicopter_dag_2_B.Gain1[0] = helicopter_dag_2_P.Gain1_Gain *
    helicopter_dag_2_B.Gain;
  helicopter_dag_2_B.Gain1[1] = helicopter_dag_2_P.Gain1_Gain *
    helicopter_dag_2_B.Gain_d;
  helicopter_dag_2_B.Gain1[2] = helicopter_dag_2_P.Gain1_Gain *
    helicopter_dag_2_B.Gain_i;
  helicopter_dag_2_B.Gain1[3] = helicopter_dag_2_P.Gain1_Gain *
    helicopter_dag_2_B.Gain_b;
  helicopter_dag_2_B.Gain1[4] = helicopter_dag_2_P.Gain1_Gain *
    helicopter_dag_2_B.Sum;
  helicopter_dag_2_B.Gain1[5] = helicopter_dag_2_P.Gain1_Gain *
    helicopter_dag_2_B.Gain_dg;
  if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
    /* SignalConversion generated from: '<Root>/To File' */
    helicopter_dag_2_B.TmpSignalConversionAtToFileInpo[0] =
      helicopter_dag_2_B.FromWorkspace;
    helicopter_dag_2_B.TmpSignalConversionAtToFileInpo[1] =
      helicopter_dag_2_B.Gain1[0];

    /* ToFile: '<Root>/To File' */
    if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
      if (rtmIsMajorTimeStep(helicopter_dag_2_M) ) {
        {
          if (!(++helicopter_dag_2_DW.ToFile_IWORK.Decimation % 2) &&
              (helicopter_dag_2_DW.ToFile_IWORK.Count * (2 + 1)) + 1 < 100000000
              ) {
            FILE *fp = (FILE *) helicopter_dag_2_DW.ToFile_PWORK.FilePtr;
            if (fp != (NULL)) {
              real_T u[2 + 1];
              helicopter_dag_2_DW.ToFile_IWORK.Decimation = 0;
              u[0] = helicopter_dag_2_M->Timing.t[1];
              u[1] = helicopter_dag_2_B.TmpSignalConversionAtToFileInpo[0];
              u[2] = helicopter_dag_2_B.TmpSignalConversionAtToFileInpo[1];
              if (fwrite(u, sizeof(real_T), 2 + 1, fp) != 2 + 1) {
                rtmSetErrorStatus(helicopter_dag_2_M,
                                  "Error writing to MAT-file lambda_vs_pc.mat");
                return;
              }

              if (((++helicopter_dag_2_DW.ToFile_IWORK.Count) * (2 + 1))+1 >=
                  100000000) {
                (void)fprintf(stdout,
                              "*** The ToFile block will stop logging data before\n"
                              "    the simulation has ended, because it has reached\n"
                              "    the maximum number of elements (100000000)\n"
                              "    allowed in MAT-file lambda_vs_pc.mat.\n");
              }
            }
          }
        }
      }
    }
  }

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S6>/K_pd'
   *  Gain: '<S6>/K_pp'
   *  Sum: '<S6>/Sum2'
   *  Sum: '<S6>/Sum3'
   */
  rtb_Clock = ((helicopter_dag_2_B.FromWorkspace - helicopter_dag_2_B.Gain1[2]) *
               helicopter_dag_2_P.K_pp - helicopter_dag_2_P.K_pd *
               helicopter_dag_2_B.Gain1[3]) + helicopter_dag_2_P.Vd_ff;

  /* Integrator: '<S3>/Integrator' */
  /* Limited  Integrator  */
  if (helicopter_dag_2_X.Integrator_CSTATE >=
      helicopter_dag_2_P.Integrator_UpperSat) {
    helicopter_dag_2_X.Integrator_CSTATE =
      helicopter_dag_2_P.Integrator_UpperSat;
  } else {
    if (helicopter_dag_2_X.Integrator_CSTATE <=
        helicopter_dag_2_P.Integrator_LowerSat) {
      helicopter_dag_2_X.Integrator_CSTATE =
        helicopter_dag_2_P.Integrator_LowerSat;
    }
  }

  /* RateTransition: '<S5>/Rate Transition: y' */
  if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
    if (helicopter_dag_2_M->Timing.RateInteraction.TID1_2) {
      /* RateTransition: '<S5>/Rate Transition: y' */
      helicopter_dag_2_B.RateTransitiony =
        helicopter_dag_2_DW.RateTransitiony_Buffer0;
    }

    /* DeadZone: '<S5>/Dead Zone: y' */
    if (helicopter_dag_2_B.RateTransitiony > helicopter_dag_2_P.DeadZoney_End) {
      /* DeadZone: '<S5>/Dead Zone: x' */
      rtb_DeadZonex = helicopter_dag_2_B.RateTransitiony -
        helicopter_dag_2_P.DeadZoney_End;
    } else if (helicopter_dag_2_B.RateTransitiony >=
               helicopter_dag_2_P.DeadZoney_Start) {
      /* DeadZone: '<S5>/Dead Zone: x' */
      rtb_DeadZonex = 0.0;
    } else {
      /* DeadZone: '<S5>/Dead Zone: x' */
      rtb_DeadZonex = helicopter_dag_2_B.RateTransitiony -
        helicopter_dag_2_P.DeadZoney_Start;
    }

    /* End of DeadZone: '<S5>/Dead Zone: y' */

    /* Gain: '<S5>/Joystick_gain_y' incorporates:
     *  Gain: '<S5>/Gain: y'
     */
    helicopter_dag_2_B.Joystick_gain_y = helicopter_dag_2_P.Gainy_Gain *
      rtb_DeadZonex * helicopter_dag_2_P.Joystick_gain_y_Gain;
  }

  /* End of RateTransition: '<S5>/Rate Transition: y' */

  /* Sum: '<S3>/Sum' */
  rtb_Derivative = helicopter_dag_2_B.Joystick_gain_y -
    helicopter_dag_2_B.Gain1[4];

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Integrator: '<S3>/Integrator'
   *  Sum: '<S3>/Sum1'
   */
  rtb_Backgain = ((helicopter_dag_2_P.K_ep * rtb_Derivative +
                   helicopter_dag_2_X.Integrator_CSTATE) -
                  helicopter_dag_2_P.K_ed * helicopter_dag_2_B.Gain1[5]) +
    helicopter_dag_2_P.Vs_ff;

  /* If: '<S3>/If' incorporates:
   *  Clock: '<S3>/Clock'
   *  Gain: '<S3>/K_ei'
   *  Inport: '<S7>/In1'
   */
  if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
    rtAction = (int8_T)!(helicopter_dag_2_M->Timing.t[0] >= 2.0);
    helicopter_dag_2_DW.If_ActiveSubsystem = rtAction;
  } else {
    rtAction = helicopter_dag_2_DW.If_ActiveSubsystem;
  }

  if (rtAction == 0) {
    /* Outputs for IfAction SubSystem: '<S3>/If Action Subsystem' incorporates:
     *  ActionPort: '<S7>/Action Port'
     */
    helicopter_dag_2_B.In1 = helicopter_dag_2_P.K_ei * rtb_Derivative;
    if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
      srUpdateBC(helicopter_dag_2_DW.IfActionSubsystem_SubsysRanBC);
    }

    /* End of Outputs for SubSystem: '<S3>/If Action Subsystem' */
  }

  /* End of If: '<S3>/If' */
  if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  rtb_Derivative = helicopter_dag_2_M->Timing.t[0];
  if ((helicopter_dag_2_DW.TimeStampA >= rtb_Derivative) &&
      (helicopter_dag_2_DW.TimeStampB >= rtb_Derivative)) {
    rtb_Derivative = 0.0;
  } else {
    lastTime = helicopter_dag_2_DW.TimeStampA;
    lastU = &helicopter_dag_2_DW.LastUAtTimeA;
    if (helicopter_dag_2_DW.TimeStampA < helicopter_dag_2_DW.TimeStampB) {
      if (helicopter_dag_2_DW.TimeStampB < rtb_Derivative) {
        lastTime = helicopter_dag_2_DW.TimeStampB;
        lastU = &helicopter_dag_2_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter_dag_2_DW.TimeStampA >= rtb_Derivative) {
        lastTime = helicopter_dag_2_DW.TimeStampB;
        lastU = &helicopter_dag_2_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter_dag_2_B.PitchCounttorad - *lastU) /
      (rtb_Derivative - lastTime);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S12>/Gain' */
  helicopter_dag_2_B.Gain_l = helicopter_dag_2_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Derivative = (rtb_Backgain - rtb_Clock) * helicopter_dag_2_P.Backgain_Gain;

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Derivative > helicopter_dag_2_P.BackmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Back motor: Saturation' */
    helicopter_dag_2_B.BackmotorSaturation =
      helicopter_dag_2_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Derivative < helicopter_dag_2_P.BackmotorSaturation_LowerSat) {
    /* Saturate: '<S4>/Back motor: Saturation' */
    helicopter_dag_2_B.BackmotorSaturation =
      helicopter_dag_2_P.BackmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Back motor: Saturation' */
    helicopter_dag_2_B.BackmotorSaturation = rtb_Derivative;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Derivative = (rtb_Clock + rtb_Backgain) *
    helicopter_dag_2_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Derivative > helicopter_dag_2_P.FrontmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Front motor: Saturation' */
    helicopter_dag_2_B.FrontmotorSaturation =
      helicopter_dag_2_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Derivative < helicopter_dag_2_P.FrontmotorSaturation_LowerSat)
  {
    /* Saturate: '<S4>/Front motor: Saturation' */
    helicopter_dag_2_B.FrontmotorSaturation =
      helicopter_dag_2_P.FrontmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Front motor: Saturation' */
    helicopter_dag_2_B.FrontmotorSaturation = rtb_Derivative;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helicopter_dag_2/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter_dag_2_DW.HILWriteAnalog_Buffer[0] =
        helicopter_dag_2_B.FrontmotorSaturation;
      helicopter_dag_2_DW.HILWriteAnalog_Buffer[1] =
        helicopter_dag_2_B.BackmotorSaturation;
      result = hil_write_analog(helicopter_dag_2_DW.HILInitialize_Card,
        helicopter_dag_2_P.HILWriteAnalog_channels, 2,
        &helicopter_dag_2_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
      }
    }

    /* RateTransition: '<S5>/Rate Transition: x' */
    if (helicopter_dag_2_M->Timing.RateInteraction.TID1_2) {
      /* RateTransition: '<S5>/Rate Transition: x' */
      helicopter_dag_2_B.RateTransitionx =
        helicopter_dag_2_DW.RateTransitionx_Buffer0;
    }

    /* End of RateTransition: '<S5>/Rate Transition: x' */

    /* DeadZone: '<S5>/Dead Zone: x' */
    if (helicopter_dag_2_B.RateTransitionx > helicopter_dag_2_P.DeadZonex_End) {
      /* DeadZone: '<S5>/Dead Zone: x' */
      rtb_DeadZonex = helicopter_dag_2_B.RateTransitionx -
        helicopter_dag_2_P.DeadZonex_End;
    } else if (helicopter_dag_2_B.RateTransitionx >=
               helicopter_dag_2_P.DeadZonex_Start) {
      /* DeadZone: '<S5>/Dead Zone: x' */
      rtb_DeadZonex = 0.0;
    } else {
      /* DeadZone: '<S5>/Dead Zone: x' */
      rtb_DeadZonex = helicopter_dag_2_B.RateTransitionx -
        helicopter_dag_2_P.DeadZonex_Start;
    }

    /* End of DeadZone: '<S5>/Dead Zone: x' */

    /* Gain: '<S5>/Joystick_gain_x' incorporates:
     *  Gain: '<S5>/Gain: x'
     */
    helicopter_dag_2_B.Joystick_gain_x = helicopter_dag_2_P.Gainx_Gain *
      rtb_DeadZonex * helicopter_dag_2_P.Joystick_gain_x_Gain;
  }
}

/* Model update function for TID0 */
void helicopter_dag_2_update0(void)    /* Sample time: [0.0s, 0.0s] */
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (helicopter_dag_2_DW.TimeStampA == (rtInf)) {
    helicopter_dag_2_DW.TimeStampA = helicopter_dag_2_M->Timing.t[0];
    lastU = &helicopter_dag_2_DW.LastUAtTimeA;
  } else if (helicopter_dag_2_DW.TimeStampB == (rtInf)) {
    helicopter_dag_2_DW.TimeStampB = helicopter_dag_2_M->Timing.t[0];
    lastU = &helicopter_dag_2_DW.LastUAtTimeB;
  } else if (helicopter_dag_2_DW.TimeStampA < helicopter_dag_2_DW.TimeStampB) {
    helicopter_dag_2_DW.TimeStampA = helicopter_dag_2_M->Timing.t[0];
    lastU = &helicopter_dag_2_DW.LastUAtTimeA;
  } else {
    helicopter_dag_2_DW.TimeStampB = helicopter_dag_2_M->Timing.t[0];
    lastU = &helicopter_dag_2_DW.LastUAtTimeB;
  }

  *lastU = helicopter_dag_2_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter_dag_2_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter_dag_2_M->solverInfo);
  }

  /* Update absolute time */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helicopter_dag_2_M->Timing.clockTick0)) {
    ++helicopter_dag_2_M->Timing.clockTickH0;
  }

  helicopter_dag_2_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicopter_dag_2_M->solverInfo);

  /* Update absolute time */
  /* The "clockTick1" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick1"
   * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick1 and the high bits
   * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helicopter_dag_2_M->Timing.clockTick1)) {
    ++helicopter_dag_2_M->Timing.clockTickH1;
  }

  helicopter_dag_2_M->Timing.t[1] = helicopter_dag_2_M->Timing.clockTick1 *
    helicopter_dag_2_M->Timing.stepSize1 +
    helicopter_dag_2_M->Timing.clockTickH1 *
    helicopter_dag_2_M->Timing.stepSize1 * 4294967296.0;
}

/* Derivatives for root system: '<Root>' */
void helicopter_dag_2_derivatives(void)
{
  XDot_helicopter_dag_2_T *_rtXdot;
  boolean_T lsat;
  boolean_T usat;
  _rtXdot = ((XDot_helicopter_dag_2_T *) helicopter_dag_2_M->derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter_dag_2_P.TravelTransferFcn_A *
    helicopter_dag_2_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter_dag_2_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_dag_2_P.PitchTransferFcn_A *
    helicopter_dag_2_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_dag_2_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter_dag_2_P.ElevationTransferFcn_A *
    helicopter_dag_2_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter_dag_2_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  lsat = (helicopter_dag_2_X.Integrator_CSTATE <=
          helicopter_dag_2_P.Integrator_LowerSat);
  usat = (helicopter_dag_2_X.Integrator_CSTATE >=
          helicopter_dag_2_P.Integrator_UpperSat);
  if (((!lsat) && (!usat)) || (lsat && (helicopter_dag_2_B.In1 > 0.0)) || (usat &&
       (helicopter_dag_2_B.In1 < 0.0))) {
    _rtXdot->Integrator_CSTATE = helicopter_dag_2_B.In1;
  } else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S3>/Integrator' */
}

/* Model output function for TID2 */
void helicopter_dag_2_output2(void)    /* Sample time: [0.01s, 0.0s] */
{
  /* local block i/o variables */
  real_T rtb_GameController_o4;
  real_T rtb_GameController_o5;

  /* S-Function (game_controller_block): '<S5>/Game Controller' */

  /* S-Function Block: helicopter_dag_2/Joystick/Game Controller (game_controller_block) */
  {
    if (helicopter_dag_2_P.GameController_Enabled) {
      t_game_controller_states state;
      t_boolean new_data;
      t_error result;
      result = game_controller_poll
        (helicopter_dag_2_DW.GameController_Controller, &state, &new_data);
      if (result < 0) {
        new_data = false;
      }

      rtb_GameController_o4 = state.x;
      rtb_GameController_o5 = state.y;
    } else {
      rtb_GameController_o4 = 0;
      rtb_GameController_o5 = 0;
    }
  }

  /* RateTransition: '<S5>/Rate Transition: x' */
  helicopter_dag_2_DW.RateTransitionx_Buffer0 = rtb_GameController_o4;

  /* RateTransition: '<S5>/Rate Transition: y' */
  helicopter_dag_2_DW.RateTransitiony_Buffer0 = rtb_GameController_o5;
}

/* Model update function for TID2 */
void helicopter_dag_2_update2(void)    /* Sample time: [0.01s, 0.0s] */
{
  /* Update absolute time */
  /* The "clockTick2" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick2"
   * and "Timing.stepSize2". Size of "clockTick2" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick2 and the high bits
   * Timing.clockTickH2. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helicopter_dag_2_M->Timing.clockTick2)) {
    ++helicopter_dag_2_M->Timing.clockTickH2;
  }

  helicopter_dag_2_M->Timing.t[2] = helicopter_dag_2_M->Timing.clockTick2 *
    helicopter_dag_2_M->Timing.stepSize2 +
    helicopter_dag_2_M->Timing.clockTickH2 *
    helicopter_dag_2_M->Timing.stepSize2 * 4294967296.0;
}

/* Model output wrapper function for compatibility with a static main program */
void helicopter_dag_2_output(int_T tid)
{
  switch (tid) {
   case 0 :
    helicopter_dag_2_output0();
    break;

   case 2 :
    helicopter_dag_2_output2();
    break;

   default :
    break;
  }
}

/* Model update wrapper function for compatibility with a static main program */
void helicopter_dag_2_update(int_T tid)
{
  switch (tid) {
   case 0 :
    helicopter_dag_2_update0();
    break;

   case 2 :
    helicopter_dag_2_update2();
    break;

   default :
    break;
  }
}

/* Model initialize function */
void helicopter_dag_2_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter_dag_2/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter_dag_2_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options
      (helicopter_dag_2_DW.HILInitialize_Card, "update_rate=normal;decimation=1",
       32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter_dag_2_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
      return;
    }

    if ((helicopter_dag_2_P.HILInitialize_AIPStart && !is_switching) ||
        (helicopter_dag_2_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter_dag_2_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = (helicopter_dag_2_P.HILInitialize_AILow);
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter_dag_2_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helicopter_dag_2_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges
        (helicopter_dag_2_DW.HILInitialize_Card,
         helicopter_dag_2_P.HILInitialize_AIChannels, 8U,
         &helicopter_dag_2_DW.HILInitialize_AIMinimums[0],
         &helicopter_dag_2_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_dag_2_P.HILInitialize_AOPStart && !is_switching) ||
        (helicopter_dag_2_P.HILInitialize_AOPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter_dag_2_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = (helicopter_dag_2_P.HILInitialize_AOLow);
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter_dag_2_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helicopter_dag_2_P.HILInitialize_AOHigh;
        }
      }

      result = hil_set_analog_output_ranges
        (helicopter_dag_2_DW.HILInitialize_Card,
         helicopter_dag_2_P.HILInitialize_AOChannels, 8U,
         &helicopter_dag_2_DW.HILInitialize_AOMinimums[0],
         &helicopter_dag_2_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_dag_2_P.HILInitialize_AOStart && !is_switching) ||
        (helicopter_dag_2_P.HILInitialize_AOEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_dag_2_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_dag_2_P.HILInitialize_AOInitial;
        }
      }

      result = hil_write_analog(helicopter_dag_2_DW.HILInitialize_Card,
        helicopter_dag_2_P.HILInitialize_AOChannels, 8U,
        &helicopter_dag_2_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_dag_2_P.HILInitialize_AOReset) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_dag_2_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_dag_2_P.HILInitialize_AOWatchdog;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter_dag_2_DW.HILInitialize_Card,
         helicopter_dag_2_P.HILInitialize_AOChannels, 8U,
         &helicopter_dag_2_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_dag_2_P.HILInitialize_EIPStart && !is_switching) ||
        (helicopter_dag_2_P.HILInitialize_EIPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter_dag_2_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter_dag_2_P.HILInitialize_EIQuadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicopter_dag_2_DW.HILInitialize_Card,
         helicopter_dag_2_P.HILInitialize_EIChannels, 8U,
         (t_encoder_quadrature_mode *)
         &helicopter_dag_2_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_dag_2_P.HILInitialize_EIStart && !is_switching) ||
        (helicopter_dag_2_P.HILInitialize_EIEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter_dag_2_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = helicopter_dag_2_P.HILInitialize_EIInitial;
        }
      }

      result = hil_set_encoder_counts(helicopter_dag_2_DW.HILInitialize_Card,
        helicopter_dag_2_P.HILInitialize_EIChannels, 8U,
        &helicopter_dag_2_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_dag_2_P.HILInitialize_POPStart && !is_switching) ||
        (helicopter_dag_2_P.HILInitialize_POPEnter && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter_dag_2_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter_dag_2_P.HILInitialize_POModes;
        }
      }

      result = hil_set_pwm_mode(helicopter_dag_2_DW.HILInitialize_Card,
        helicopter_dag_2_P.HILInitialize_POChannels, 8U, (t_pwm_mode *)
        &helicopter_dag_2_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_POChannels =
          helicopter_dag_2_P.HILInitialize_POChannels;
        int32_T *dw_POModeValues =
          &helicopter_dag_2_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE ||
              dw_POModeValues[i1] == PWM_RAW_MODE) {
            helicopter_dag_2_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = (p_HILInitialize_POChannels[i1]);
            helicopter_dag_2_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = helicopter_dag_2_P.HILInitialize_POFrequency;
            num_duty_cycle_modes++;
          } else {
            helicopter_dag_2_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = (p_HILInitialize_POChannels[i1]);
            helicopter_dag_2_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] =
              helicopter_dag_2_P.HILInitialize_POFrequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter_dag_2_DW.HILInitialize_Card,
          &helicopter_dag_2_DW.HILInitialize_POSortedChans[0],
          num_duty_cycle_modes,
          &helicopter_dag_2_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter_dag_2_DW.HILInitialize_Card,
          &helicopter_dag_2_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter_dag_2_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter_dag_2_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter_dag_2_P.HILInitialize_POConfiguration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicopter_dag_2_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicopter_dag_2_P.HILInitialize_POAlignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter_dag_2_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicopter_dag_2_P.HILInitialize_POPolarity;
        }
      }

      result = hil_set_pwm_configuration(helicopter_dag_2_DW.HILInitialize_Card,
        helicopter_dag_2_P.HILInitialize_POChannels, 8U,
        (t_pwm_configuration *) &helicopter_dag_2_DW.HILInitialize_POModeValues
        [0],
        (t_pwm_alignment *) &helicopter_dag_2_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicopter_dag_2_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &helicopter_dag_2_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = helicopter_dag_2_P.HILInitialize_POLeading;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter_dag_2_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_dag_2_P.HILInitialize_POTrailing;
        }
      }

      result = hil_set_pwm_deadband(helicopter_dag_2_DW.HILInitialize_Card,
        helicopter_dag_2_P.HILInitialize_POChannels, 8U,
        &helicopter_dag_2_DW.HILInitialize_POSortedFreqs[0],
        &helicopter_dag_2_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_dag_2_P.HILInitialize_POStart && !is_switching) ||
        (helicopter_dag_2_P.HILInitialize_POEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_dag_2_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_dag_2_P.HILInitialize_POInitial;
        }
      }

      result = hil_write_pwm(helicopter_dag_2_DW.HILInitialize_Card,
        helicopter_dag_2_P.HILInitialize_POChannels, 8U,
        &helicopter_dag_2_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_dag_2_P.HILInitialize_POReset) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_dag_2_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_dag_2_P.HILInitialize_POWatchdog;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter_dag_2_DW.HILInitialize_Card,
         helicopter_dag_2_P.HILInitialize_POChannels, 8U,
         &helicopter_dag_2_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter_dag_2/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (helicopter_dag_2_DW.HILInitialize_Card,
       helicopter_dag_2_P.HILReadEncoderTimebase_SamplesI,
       helicopter_dag_2_P.HILReadEncoderTimebase_Channels, 3,
       &helicopter_dag_2_DW.HILReadEncoderTimebase_Task);
    if (result >= 0) {
      result = hil_task_set_buffer_overflow_mode
        (helicopter_dag_2_DW.HILReadEncoderTimebase_Task,
         (t_buffer_overflow_mode)
         (helicopter_dag_2_P.HILReadEncoderTimebase_Overflow - 1));
    }

    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
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

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559828927,
      0.52359877559828771, 0.52359877559828572, 0.52359877559828327,
      0.52359877559827972, 0.52359877559827528, 0.52359877559826928,
      0.52359877559826062, 0.52359877559824752, 0.52359877559822976,
      0.52359877559816226, 0.52359877559762713, 0.52359877559764423,
      0.0086840231302012239, -0.52359877559751444, -0.52359877559792189,
      -0.52359877559801293, -0.52359877559804535, -0.52359877559805368,
      -0.52359877559804546, -0.5235987755980227, -0.52359877559797,
      -0.523598775597857, -0.52359877559763535, -0.52359877559690049,
      -0.52359877559658374, -0.52359877555341372, -0.52359877556069856,
      -0.47307630323645, -0.19354597354305425, 0.002217996219539109,
      0.1260562620884359, 0.19193445751144478, 0.21420857674488203,
      0.20629374383124066, 0.17978711806924108, 0.14401542121825328,
      0.10593114545041349, 0.070262864966350258, 0.039825795824210086,
      0.015911104794662623, -0.0013096278980783138, -0.012409256122354329,
      -0.018386086855065265, -0.020413359368047113, -0.019658509196256557,
      -0.017166135556775486, -0.013794311282383731, -0.01019209914157615,
      -0.0068061941804854476, -0.0039058654239925161, -0.0016172829337297578,
      3.9533089792431753E-5, 0.0011158157986229655, 0.0017040458391527347,
      0.0019141884107221063, 0.0018564193874759383, 0.0016298405148635897,
      0.0013162481380808444, 0.000977831161999343, 0.00065766529427901954,
      0.00038198028433289455, 0.00016335250955490377, 4.1760433170390243E-6,
      -0.00010003791567236853, -0.00015780928249831128, -0.00017941942045685,
      -0.00017525198741918935, -0.00015470069632184114, -0.00012555787225554393,
      -9.3778727479842949E-5, -6.3514605156633763E-5, -3.7318216213111377E-5,
      -1.644013488000251E-5, -1.1546240067383806E-6, 8.92857465106811E-6,
      1.4593332279710935E-5, 1.6800451453247156E-5, 1.6528392787762414E-5,
      1.4667640291565753E-5, 1.1960786202913631E-5, 8.9785153447596855E-6,
      6.12142823974704E-6, 3.6385087023305118E-6, 1.6545425692848781E-6,
      2.0055718130684852E-7, -7.5689293232983346E-7, -1.2901942513243725E-6,
      -1.4901980521786484E-6, -1.4511124336591408E-6, -1.2603855124027064E-6,
      -9.92775771146981E-7, -7.0759256121100123E-7, -4.4804221477345862E-7,
      -2.4169700441056818E-7, -1.0130744676484227E-7, -2.5550311200106535E-8,
      0.0, 0.0, 1.1102230246251565E-16, 1.1102230246251565E-16, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0 } ;

    helicopter_dag_2_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter_dag_2_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helicopter_dag_2_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "lambda_vs_pc.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helicopter_dag_2_M,
                        "Error creating .mat file lambda_vs_pc.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp, 2 + 1, 0, "ans")) {
      rtmSetErrorStatus(helicopter_dag_2_M,
                        "Error writing mat file header to file lambda_vs_pc.mat");
      return;
    }

    helicopter_dag_2_DW.ToFile_IWORK.Count = 0;
    helicopter_dag_2_DW.ToFile_IWORK.Decimation = -1;
    helicopter_dag_2_DW.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for RateTransition: '<S5>/Rate Transition: y' */
  helicopter_dag_2_B.RateTransitiony =
    helicopter_dag_2_P.RateTransitiony_InitialConditio;

  /* Start for If: '<S3>/If' */
  helicopter_dag_2_DW.If_ActiveSubsystem = -1;

  /* Start for RateTransition: '<S5>/Rate Transition: x' */
  helicopter_dag_2_B.RateTransitionx =
    helicopter_dag_2_P.RateTransitionx_InitialConditio;

  /* Start for S-Function (game_controller_block): '<S5>/Game Controller' */

  /* S-Function Block: helicopter_dag_2/Joystick/Game Controller (game_controller_block) */
  {
    if (helicopter_dag_2_P.GameController_Enabled) {
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
        (helicopter_dag_2_P.GameController_ControllerNumber,
         helicopter_dag_2_P.GameController_BufferSize, deadzone, saturation,
         helicopter_dag_2_P.GameController_AutoCenter, 0, 1.0,
         &helicopter_dag_2_DW.GameController_Controller);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
      }
    }
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helicopter_dag_2_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helicopter_dag_2_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helicopter_dag_2_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  helicopter_dag_2_X.Integrator_CSTATE = helicopter_dag_2_P.Integrator_IC;

  /* InitializeConditions for RateTransition: '<S5>/Rate Transition: y' */
  helicopter_dag_2_DW.RateTransitiony_Buffer0 =
    helicopter_dag_2_P.RateTransitiony_InitialConditio;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  helicopter_dag_2_DW.TimeStampA = (rtInf);
  helicopter_dag_2_DW.TimeStampB = (rtInf);

  /* InitializeConditions for RateTransition: '<S5>/Rate Transition: x' */
  helicopter_dag_2_DW.RateTransitionx_Buffer0 =
    helicopter_dag_2_P.RateTransitionx_InitialConditio;
}

/* Model terminate function */
void helicopter_dag_2_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter_dag_2/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter_dag_2_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter_dag_2_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter_dag_2_P.HILInitialize_AOTerminate && !is_switching) ||
        (helicopter_dag_2_P.HILInitialize_AOExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_dag_2_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_dag_2_P.HILInitialize_AOFinal;
        }
      }

      num_final_analog_outputs = 8U;
    } else {
      num_final_analog_outputs = 0;
    }

    if ((helicopter_dag_2_P.HILInitialize_POTerminate && !is_switching) ||
        (helicopter_dag_2_P.HILInitialize_POExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_dag_2_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_dag_2_P.HILInitialize_POFinal;
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
      result = hil_write(helicopter_dag_2_DW.HILInitialize_Card
                         , helicopter_dag_2_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , helicopter_dag_2_P.HILInitialize_POChannels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter_dag_2_DW.HILInitialize_AOVoltages[0]
                         , &helicopter_dag_2_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helicopter_dag_2_DW.HILInitialize_Card,
            helicopter_dag_2_P.HILInitialize_AOChannels,
            num_final_analog_outputs,
            &helicopter_dag_2_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter_dag_2_DW.HILInitialize_Card,
            helicopter_dag_2_P.HILInitialize_POChannels, num_final_pwm_outputs,
            &helicopter_dag_2_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_dag_2_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter_dag_2_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter_dag_2_DW.HILInitialize_Card);
    hil_close(helicopter_dag_2_DW.HILInitialize_Card);
    helicopter_dag_2_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helicopter_dag_2_DW.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "lambda_vs_pc.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_dag_2_M,
                          "Error closing MAT-file lambda_vs_pc.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helicopter_dag_2_M,
                          "Error reopening MAT-file lambda_vs_pc.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2 + 1,
           helicopter_dag_2_DW.ToFile_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helicopter_dag_2_M,
                          "Error writing header for ans to MAT-file lambda_vs_pc.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_dag_2_M,
                          "Error closing MAT-file lambda_vs_pc.mat");
        return;
      }

      helicopter_dag_2_DW.ToFile_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for S-Function (game_controller_block): '<S5>/Game Controller' */

  /* S-Function Block: helicopter_dag_2/Joystick/Game Controller (game_controller_block) */
  {
    if (helicopter_dag_2_P.GameController_Enabled) {
      game_controller_close(helicopter_dag_2_DW.GameController_Controller);
      helicopter_dag_2_DW.GameController_Controller = NULL;
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
  if (tid == 1)
    tid = 0;
  helicopter_dag_2_output(tid);
}

void MdlUpdate(int_T tid)
{
  if (tid == 1)
    tid = 0;
  helicopter_dag_2_update(tid);
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
  helicopter_dag_2_initialize();
}

void MdlTerminate(void)
{
  helicopter_dag_2_terminate();
}

/* Registration function */
RT_MODEL_helicopter_dag_2_T *helicopter_dag_2(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter_dag_2_P.Integrator_UpperSat = rtInf;
  helicopter_dag_2_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter_dag_2_M, 0,
                sizeof(RT_MODEL_helicopter_dag_2_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter_dag_2_M->solverInfo,
                          &helicopter_dag_2_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter_dag_2_M->solverInfo, &rtmGetTPtr(helicopter_dag_2_M));
    rtsiSetStepSizePtr(&helicopter_dag_2_M->solverInfo,
                       &helicopter_dag_2_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter_dag_2_M->solverInfo, &helicopter_dag_2_M->derivs);
    rtsiSetContStatesPtr(&helicopter_dag_2_M->solverInfo, (real_T **)
                         &helicopter_dag_2_M->contStates);
    rtsiSetNumContStatesPtr(&helicopter_dag_2_M->solverInfo,
      &helicopter_dag_2_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&helicopter_dag_2_M->solverInfo,
      &helicopter_dag_2_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&helicopter_dag_2_M->solverInfo,
      &helicopter_dag_2_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&helicopter_dag_2_M->solverInfo,
      &helicopter_dag_2_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&helicopter_dag_2_M->solverInfo, (&rtmGetErrorStatus
      (helicopter_dag_2_M)));
    rtsiSetRTModelPtr(&helicopter_dag_2_M->solverInfo, helicopter_dag_2_M);
  }

  rtsiSetSimTimeStep(&helicopter_dag_2_M->solverInfo, MAJOR_TIME_STEP);
  helicopter_dag_2_M->intgData.f[0] = helicopter_dag_2_M->odeF[0];
  helicopter_dag_2_M->contStates = ((real_T *) &helicopter_dag_2_X);
  rtsiSetSolverData(&helicopter_dag_2_M->solverInfo, (void *)
                    &helicopter_dag_2_M->intgData);
  rtsiSetSolverName(&helicopter_dag_2_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter_dag_2_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    mdlTsMap[2] = 2;
    helicopter_dag_2_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter_dag_2_M->Timing.sampleTimes =
      (&helicopter_dag_2_M->Timing.sampleTimesArray[0]);
    helicopter_dag_2_M->Timing.offsetTimes =
      (&helicopter_dag_2_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter_dag_2_M->Timing.sampleTimes[0] = (0.0);
    helicopter_dag_2_M->Timing.sampleTimes[1] = (0.002);
    helicopter_dag_2_M->Timing.sampleTimes[2] = (0.01);

    /* task offsets */
    helicopter_dag_2_M->Timing.offsetTimes[0] = (0.0);
    helicopter_dag_2_M->Timing.offsetTimes[1] = (0.0);
    helicopter_dag_2_M->Timing.offsetTimes[2] = (0.0);
  }

  rtmSetTPtr(helicopter_dag_2_M, &helicopter_dag_2_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter_dag_2_M->Timing.sampleHitArray;
    int_T *mdlPerTaskSampleHits =
      helicopter_dag_2_M->Timing.perTaskSampleHitsArray;
    helicopter_dag_2_M->Timing.perTaskSampleHits = (&mdlPerTaskSampleHits[0]);
    mdlSampleHits[0] = 1;
    helicopter_dag_2_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter_dag_2_M, 20.0);
  helicopter_dag_2_M->Timing.stepSize0 = 0.002;
  helicopter_dag_2_M->Timing.stepSize1 = 0.002;
  helicopter_dag_2_M->Timing.stepSize2 = 0.01;

  /* External mode info */
  helicopter_dag_2_M->Sizes.checksums[0] = (266776760U);
  helicopter_dag_2_M->Sizes.checksums[1] = (680672581U);
  helicopter_dag_2_M->Sizes.checksums[2] = (2076109412U);
  helicopter_dag_2_M->Sizes.checksums[3] = (2763510147U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[2];
    helicopter_dag_2_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = (sysRanDType *)
      &helicopter_dag_2_DW.IfActionSubsystem_SubsysRanBC;
    rteiSetModelMappingInfoPtr(helicopter_dag_2_M->extModeInfo,
      &helicopter_dag_2_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter_dag_2_M->extModeInfo,
                        helicopter_dag_2_M->Sizes.checksums);
    rteiSetTPtr(helicopter_dag_2_M->extModeInfo, rtmGetTPtr(helicopter_dag_2_M));
  }

  helicopter_dag_2_M->solverInfoPtr = (&helicopter_dag_2_M->solverInfo);
  helicopter_dag_2_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter_dag_2_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter_dag_2_M->solverInfo, SOLVER_MODE_MULTITASKING);

  /* block I/O */
  helicopter_dag_2_M->blockIO = ((void *) &helicopter_dag_2_B);

  {
    int32_T i;
    for (i = 0; i < 6; i++) {
      helicopter_dag_2_B.Gain1[i] = 0.0;
    }

    helicopter_dag_2_B.FromWorkspace = 0.0;
    helicopter_dag_2_B.TravelCounttorad = 0.0;
    helicopter_dag_2_B.Gain = 0.0;
    helicopter_dag_2_B.Gain_d = 0.0;
    helicopter_dag_2_B.PitchCounttorad = 0.0;
    helicopter_dag_2_B.Gain_i = 0.0;
    helicopter_dag_2_B.Gain_b = 0.0;
    helicopter_dag_2_B.ElevationCounttorad = 0.0;
    helicopter_dag_2_B.Gain_e = 0.0;
    helicopter_dag_2_B.Sum = 0.0;
    helicopter_dag_2_B.Gain_dg = 0.0;
    helicopter_dag_2_B.TmpSignalConversionAtToFileInpo[0] = 0.0;
    helicopter_dag_2_B.TmpSignalConversionAtToFileInpo[1] = 0.0;
    helicopter_dag_2_B.RateTransitiony = 0.0;
    helicopter_dag_2_B.Joystick_gain_y = 0.0;
    helicopter_dag_2_B.Gain_l = 0.0;
    helicopter_dag_2_B.BackmotorSaturation = 0.0;
    helicopter_dag_2_B.FrontmotorSaturation = 0.0;
    helicopter_dag_2_B.RateTransitionx = 0.0;
    helicopter_dag_2_B.Joystick_gain_x = 0.0;
    helicopter_dag_2_B.In1 = 0.0;
  }

  /* parameters */
  helicopter_dag_2_M->defaultParam = ((real_T *)&helicopter_dag_2_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter_dag_2_X;
    helicopter_dag_2_M->contStates = (x);
    (void) memset((void *)&helicopter_dag_2_X, 0,
                  sizeof(X_helicopter_dag_2_T));
  }

  /* states (dwork) */
  helicopter_dag_2_M->dwork = ((void *) &helicopter_dag_2_DW);
  (void) memset((void *)&helicopter_dag_2_DW, 0,
                sizeof(DW_helicopter_dag_2_T));

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_2_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_2_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_2_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_2_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_2_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_2_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_2_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_dag_2_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter_dag_2_DW.RateTransitiony_Buffer0 = 0.0;
  helicopter_dag_2_DW.TimeStampA = 0.0;
  helicopter_dag_2_DW.LastUAtTimeA = 0.0;
  helicopter_dag_2_DW.TimeStampB = 0.0;
  helicopter_dag_2_DW.LastUAtTimeB = 0.0;
  helicopter_dag_2_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter_dag_2_DW.HILWriteAnalog_Buffer[1] = 0.0;
  helicopter_dag_2_DW.RateTransitionx_Buffer0 = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter_dag_2_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 17;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter_dag_2_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter_dag_2_M->Sizes.numPeriodicContStates = (0);
                                      /* Number of periodic continuous states */
  helicopter_dag_2_M->Sizes.numY = (0);/* Number of model outputs */
  helicopter_dag_2_M->Sizes.numU = (0);/* Number of model inputs */
  helicopter_dag_2_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter_dag_2_M->Sizes.numSampTimes = (3);/* Number of sample times */
  helicopter_dag_2_M->Sizes.numBlocks = (70);/* Number of blocks */
  helicopter_dag_2_M->Sizes.numBlockIO = (21);/* Number of block outputs */
  helicopter_dag_2_M->Sizes.numBlockPrms = (157);/* Sum of parameter "widths" */
  return helicopter_dag_2_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
