  function targMap = targDataMap(),

  ;%***********************
  ;% Create Parameter Map *
  ;%***********************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 8;
    sectIdxOffset = 0;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc paramMap
    ;%
    paramMap.nSections           = nTotSects;
    paramMap.sectIdxOffset       = sectIdxOffset;
      paramMap.sections(nTotSects) = dumSection; %prealloc
    paramMap.nTotData            = -1;
    
    ;%
    ;% Auto data (helicopter_dag_2_P)
    ;%
      section.nData     = 9;
      section.data(9)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_P.K_ed
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% helicopter_dag_2_P.K_ei
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% helicopter_dag_2_P.K_ep
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 2;
	
	  ;% helicopter_dag_2_P.K_pd
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 3;
	
	  ;% helicopter_dag_2_P.K_pp
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 4;
	
	  ;% helicopter_dag_2_P.Vd_ff
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 5;
	
	  ;% helicopter_dag_2_P.Vs_ff
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 6;
	
	  ;% helicopter_dag_2_P.elevation_gain
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 7;
	
	  ;% helicopter_dag_2_P.travel_gain
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 8;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(1) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_P.HILWriteAnalog_channels
	  section.data(1).logicalSrcIdx = 9;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(2) = section;
      clear section
      
      section.nData     = 57;
      section.data(57)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_P.HILInitialize_OOTerminate
	  section.data(1).logicalSrcIdx = 10;
	  section.data(1).dtTransOffset = 0;
	
	  ;% helicopter_dag_2_P.HILInitialize_OOExit
	  section.data(2).logicalSrcIdx = 11;
	  section.data(2).dtTransOffset = 1;
	
	  ;% helicopter_dag_2_P.HILInitialize_OOStart
	  section.data(3).logicalSrcIdx = 12;
	  section.data(3).dtTransOffset = 2;
	
	  ;% helicopter_dag_2_P.HILInitialize_OOEnter
	  section.data(4).logicalSrcIdx = 13;
	  section.data(4).dtTransOffset = 3;
	
	  ;% helicopter_dag_2_P.HILInitialize_AOFinal
	  section.data(5).logicalSrcIdx = 14;
	  section.data(5).dtTransOffset = 4;
	
	  ;% helicopter_dag_2_P.HILInitialize_POFinal
	  section.data(6).logicalSrcIdx = 15;
	  section.data(6).dtTransOffset = 5;
	
	  ;% helicopter_dag_2_P.HILInitialize_AIHigh
	  section.data(7).logicalSrcIdx = 16;
	  section.data(7).dtTransOffset = 6;
	
	  ;% helicopter_dag_2_P.HILInitialize_AILow
	  section.data(8).logicalSrcIdx = 17;
	  section.data(8).dtTransOffset = 7;
	
	  ;% helicopter_dag_2_P.HILInitialize_AOHigh
	  section.data(9).logicalSrcIdx = 18;
	  section.data(9).dtTransOffset = 8;
	
	  ;% helicopter_dag_2_P.HILInitialize_AOLow
	  section.data(10).logicalSrcIdx = 19;
	  section.data(10).dtTransOffset = 9;
	
	  ;% helicopter_dag_2_P.HILInitialize_AOInitial
	  section.data(11).logicalSrcIdx = 20;
	  section.data(11).dtTransOffset = 10;
	
	  ;% helicopter_dag_2_P.HILInitialize_AOWatchdog
	  section.data(12).logicalSrcIdx = 21;
	  section.data(12).dtTransOffset = 11;
	
	  ;% helicopter_dag_2_P.HILInitialize_POFrequency
	  section.data(13).logicalSrcIdx = 22;
	  section.data(13).dtTransOffset = 12;
	
	  ;% helicopter_dag_2_P.HILInitialize_POLeading
	  section.data(14).logicalSrcIdx = 23;
	  section.data(14).dtTransOffset = 13;
	
	  ;% helicopter_dag_2_P.HILInitialize_POTrailing
	  section.data(15).logicalSrcIdx = 24;
	  section.data(15).dtTransOffset = 14;
	
	  ;% helicopter_dag_2_P.HILInitialize_POInitial
	  section.data(16).logicalSrcIdx = 25;
	  section.data(16).dtTransOffset = 15;
	
	  ;% helicopter_dag_2_P.HILInitialize_POWatchdog
	  section.data(17).logicalSrcIdx = 26;
	  section.data(17).dtTransOffset = 16;
	
	  ;% helicopter_dag_2_P.TravelCounttorad_Gain
	  section.data(18).logicalSrcIdx = 27;
	  section.data(18).dtTransOffset = 17;
	
	  ;% helicopter_dag_2_P.Gain_Gain
	  section.data(19).logicalSrcIdx = 28;
	  section.data(19).dtTransOffset = 18;
	
	  ;% helicopter_dag_2_P.TravelTransferFcn_A
	  section.data(20).logicalSrcIdx = 29;
	  section.data(20).dtTransOffset = 19;
	
	  ;% helicopter_dag_2_P.TravelTransferFcn_C
	  section.data(21).logicalSrcIdx = 30;
	  section.data(21).dtTransOffset = 20;
	
	  ;% helicopter_dag_2_P.TravelTransferFcn_D
	  section.data(22).logicalSrcIdx = 31;
	  section.data(22).dtTransOffset = 21;
	
	  ;% helicopter_dag_2_P.Gain_Gain_l
	  section.data(23).logicalSrcIdx = 32;
	  section.data(23).dtTransOffset = 22;
	
	  ;% helicopter_dag_2_P.PitchCounttorad_Gain
	  section.data(24).logicalSrcIdx = 33;
	  section.data(24).dtTransOffset = 23;
	
	  ;% helicopter_dag_2_P.Gain_Gain_a
	  section.data(25).logicalSrcIdx = 34;
	  section.data(25).dtTransOffset = 24;
	
	  ;% helicopter_dag_2_P.PitchTransferFcn_A
	  section.data(26).logicalSrcIdx = 35;
	  section.data(26).dtTransOffset = 25;
	
	  ;% helicopter_dag_2_P.PitchTransferFcn_C
	  section.data(27).logicalSrcIdx = 36;
	  section.data(27).dtTransOffset = 26;
	
	  ;% helicopter_dag_2_P.PitchTransferFcn_D
	  section.data(28).logicalSrcIdx = 37;
	  section.data(28).dtTransOffset = 27;
	
	  ;% helicopter_dag_2_P.Gain_Gain_ae
	  section.data(29).logicalSrcIdx = 38;
	  section.data(29).dtTransOffset = 28;
	
	  ;% helicopter_dag_2_P.ElevationCounttorad_Gain
	  section.data(30).logicalSrcIdx = 39;
	  section.data(30).dtTransOffset = 29;
	
	  ;% helicopter_dag_2_P.Gain_Gain_lv
	  section.data(31).logicalSrcIdx = 40;
	  section.data(31).dtTransOffset = 30;
	
	  ;% helicopter_dag_2_P.elavation_offsetdeg_Value
	  section.data(32).logicalSrcIdx = 41;
	  section.data(32).dtTransOffset = 31;
	
	  ;% helicopter_dag_2_P.ElevationTransferFcn_A
	  section.data(33).logicalSrcIdx = 42;
	  section.data(33).dtTransOffset = 32;
	
	  ;% helicopter_dag_2_P.ElevationTransferFcn_C
	  section.data(34).logicalSrcIdx = 43;
	  section.data(34).dtTransOffset = 33;
	
	  ;% helicopter_dag_2_P.ElevationTransferFcn_D
	  section.data(35).logicalSrcIdx = 44;
	  section.data(35).dtTransOffset = 34;
	
	  ;% helicopter_dag_2_P.Gain_Gain_n
	  section.data(36).logicalSrcIdx = 45;
	  section.data(36).dtTransOffset = 35;
	
	  ;% helicopter_dag_2_P.Gain1_Gain
	  section.data(37).logicalSrcIdx = 46;
	  section.data(37).dtTransOffset = 36;
	
	  ;% helicopter_dag_2_P.Integrator_IC
	  section.data(38).logicalSrcIdx = 47;
	  section.data(38).dtTransOffset = 37;
	
	  ;% helicopter_dag_2_P.Integrator_UpperSat
	  section.data(39).logicalSrcIdx = 48;
	  section.data(39).dtTransOffset = 38;
	
	  ;% helicopter_dag_2_P.Integrator_LowerSat
	  section.data(40).logicalSrcIdx = 49;
	  section.data(40).dtTransOffset = 39;
	
	  ;% helicopter_dag_2_P.RateTransitiony_InitialConditio
	  section.data(41).logicalSrcIdx = 50;
	  section.data(41).dtTransOffset = 40;
	
	  ;% helicopter_dag_2_P.DeadZoney_Start
	  section.data(42).logicalSrcIdx = 51;
	  section.data(42).dtTransOffset = 41;
	
	  ;% helicopter_dag_2_P.DeadZoney_End
	  section.data(43).logicalSrcIdx = 52;
	  section.data(43).dtTransOffset = 42;
	
	  ;% helicopter_dag_2_P.Gainy_Gain
	  section.data(44).logicalSrcIdx = 53;
	  section.data(44).dtTransOffset = 43;
	
	  ;% helicopter_dag_2_P.Joystick_gain_y_Gain
	  section.data(45).logicalSrcIdx = 54;
	  section.data(45).dtTransOffset = 44;
	
	  ;% helicopter_dag_2_P.Backgain_Gain
	  section.data(46).logicalSrcIdx = 55;
	  section.data(46).dtTransOffset = 45;
	
	  ;% helicopter_dag_2_P.Frontgain_Gain
	  section.data(47).logicalSrcIdx = 56;
	  section.data(47).dtTransOffset = 46;
	
	  ;% helicopter_dag_2_P.Gain_Gain_a1
	  section.data(48).logicalSrcIdx = 57;
	  section.data(48).dtTransOffset = 47;
	
	  ;% helicopter_dag_2_P.BackmotorSaturation_UpperSat
	  section.data(49).logicalSrcIdx = 58;
	  section.data(49).dtTransOffset = 48;
	
	  ;% helicopter_dag_2_P.BackmotorSaturation_LowerSat
	  section.data(50).logicalSrcIdx = 59;
	  section.data(50).dtTransOffset = 49;
	
	  ;% helicopter_dag_2_P.FrontmotorSaturation_UpperSat
	  section.data(51).logicalSrcIdx = 60;
	  section.data(51).dtTransOffset = 50;
	
	  ;% helicopter_dag_2_P.FrontmotorSaturation_LowerSat
	  section.data(52).logicalSrcIdx = 61;
	  section.data(52).dtTransOffset = 51;
	
	  ;% helicopter_dag_2_P.RateTransitionx_InitialConditio
	  section.data(53).logicalSrcIdx = 62;
	  section.data(53).dtTransOffset = 52;
	
	  ;% helicopter_dag_2_P.DeadZonex_Start
	  section.data(54).logicalSrcIdx = 63;
	  section.data(54).dtTransOffset = 53;
	
	  ;% helicopter_dag_2_P.DeadZonex_End
	  section.data(55).logicalSrcIdx = 64;
	  section.data(55).dtTransOffset = 54;
	
	  ;% helicopter_dag_2_P.Gainx_Gain
	  section.data(56).logicalSrcIdx = 65;
	  section.data(56).dtTransOffset = 55;
	
	  ;% helicopter_dag_2_P.Joystick_gain_x_Gain
	  section.data(57).logicalSrcIdx = 66;
	  section.data(57).dtTransOffset = 56;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(3) = section;
      clear section
      
      section.nData     = 8;
      section.data(8)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_P.HILInitialize_CKChannels
	  section.data(1).logicalSrcIdx = 67;
	  section.data(1).dtTransOffset = 0;
	
	  ;% helicopter_dag_2_P.HILInitialize_DOWatchdog
	  section.data(2).logicalSrcIdx = 68;
	  section.data(2).dtTransOffset = 3;
	
	  ;% helicopter_dag_2_P.HILInitialize_EIInitial
	  section.data(3).logicalSrcIdx = 69;
	  section.data(3).dtTransOffset = 4;
	
	  ;% helicopter_dag_2_P.HILInitialize_POModes
	  section.data(4).logicalSrcIdx = 70;
	  section.data(4).dtTransOffset = 5;
	
	  ;% helicopter_dag_2_P.HILInitialize_POConfiguration
	  section.data(5).logicalSrcIdx = 71;
	  section.data(5).dtTransOffset = 6;
	
	  ;% helicopter_dag_2_P.HILInitialize_POAlignment
	  section.data(6).logicalSrcIdx = 72;
	  section.data(6).dtTransOffset = 7;
	
	  ;% helicopter_dag_2_P.HILInitialize_POPolarity
	  section.data(7).logicalSrcIdx = 73;
	  section.data(7).dtTransOffset = 8;
	
	  ;% helicopter_dag_2_P.HILReadEncoderTimebase_Clock
	  section.data(8).logicalSrcIdx = 74;
	  section.data(8).dtTransOffset = 9;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(4) = section;
      clear section
      
      section.nData     = 7;
      section.data(7)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_P.HILInitialize_AIChannels
	  section.data(1).logicalSrcIdx = 75;
	  section.data(1).dtTransOffset = 0;
	
	  ;% helicopter_dag_2_P.HILInitialize_AOChannels
	  section.data(2).logicalSrcIdx = 76;
	  section.data(2).dtTransOffset = 8;
	
	  ;% helicopter_dag_2_P.HILInitialize_EIChannels
	  section.data(3).logicalSrcIdx = 77;
	  section.data(3).dtTransOffset = 16;
	
	  ;% helicopter_dag_2_P.HILInitialize_EIQuadrature
	  section.data(4).logicalSrcIdx = 78;
	  section.data(4).dtTransOffset = 24;
	
	  ;% helicopter_dag_2_P.HILInitialize_POChannels
	  section.data(5).logicalSrcIdx = 79;
	  section.data(5).dtTransOffset = 25;
	
	  ;% helicopter_dag_2_P.HILReadEncoderTimebase_Channels
	  section.data(6).logicalSrcIdx = 80;
	  section.data(6).dtTransOffset = 33;
	
	  ;% helicopter_dag_2_P.HILReadEncoderTimebase_SamplesI
	  section.data(7).logicalSrcIdx = 81;
	  section.data(7).dtTransOffset = 36;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(5) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_P.GameController_BufferSize
	  section.data(1).logicalSrcIdx = 82;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(6) = section;
      clear section
      
      section.nData     = 39;
      section.data(39)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_P.HILInitialize_Active
	  section.data(1).logicalSrcIdx = 83;
	  section.data(1).dtTransOffset = 0;
	
	  ;% helicopter_dag_2_P.HILInitialize_AOTerminate
	  section.data(2).logicalSrcIdx = 84;
	  section.data(2).dtTransOffset = 1;
	
	  ;% helicopter_dag_2_P.HILInitialize_AOExit
	  section.data(3).logicalSrcIdx = 85;
	  section.data(3).dtTransOffset = 2;
	
	  ;% helicopter_dag_2_P.HILInitialize_DOTerminate
	  section.data(4).logicalSrcIdx = 86;
	  section.data(4).dtTransOffset = 3;
	
	  ;% helicopter_dag_2_P.HILInitialize_DOExit
	  section.data(5).logicalSrcIdx = 87;
	  section.data(5).dtTransOffset = 4;
	
	  ;% helicopter_dag_2_P.HILInitialize_POTerminate
	  section.data(6).logicalSrcIdx = 88;
	  section.data(6).dtTransOffset = 5;
	
	  ;% helicopter_dag_2_P.HILInitialize_POExit
	  section.data(7).logicalSrcIdx = 89;
	  section.data(7).dtTransOffset = 6;
	
	  ;% helicopter_dag_2_P.HILInitialize_CKPStart
	  section.data(8).logicalSrcIdx = 90;
	  section.data(8).dtTransOffset = 7;
	
	  ;% helicopter_dag_2_P.HILInitialize_CKPEnter
	  section.data(9).logicalSrcIdx = 91;
	  section.data(9).dtTransOffset = 8;
	
	  ;% helicopter_dag_2_P.HILInitialize_CKStart
	  section.data(10).logicalSrcIdx = 92;
	  section.data(10).dtTransOffset = 9;
	
	  ;% helicopter_dag_2_P.HILInitialize_CKEnter
	  section.data(11).logicalSrcIdx = 93;
	  section.data(11).dtTransOffset = 10;
	
	  ;% helicopter_dag_2_P.HILInitialize_AIPStart
	  section.data(12).logicalSrcIdx = 94;
	  section.data(12).dtTransOffset = 11;
	
	  ;% helicopter_dag_2_P.HILInitialize_AIPEnter
	  section.data(13).logicalSrcIdx = 95;
	  section.data(13).dtTransOffset = 12;
	
	  ;% helicopter_dag_2_P.HILInitialize_AOPStart
	  section.data(14).logicalSrcIdx = 96;
	  section.data(14).dtTransOffset = 13;
	
	  ;% helicopter_dag_2_P.HILInitialize_AOPEnter
	  section.data(15).logicalSrcIdx = 97;
	  section.data(15).dtTransOffset = 14;
	
	  ;% helicopter_dag_2_P.HILInitialize_AOStart
	  section.data(16).logicalSrcIdx = 98;
	  section.data(16).dtTransOffset = 15;
	
	  ;% helicopter_dag_2_P.HILInitialize_AOEnter
	  section.data(17).logicalSrcIdx = 99;
	  section.data(17).dtTransOffset = 16;
	
	  ;% helicopter_dag_2_P.HILInitialize_AOReset
	  section.data(18).logicalSrcIdx = 100;
	  section.data(18).dtTransOffset = 17;
	
	  ;% helicopter_dag_2_P.HILInitialize_DOPStart
	  section.data(19).logicalSrcIdx = 101;
	  section.data(19).dtTransOffset = 18;
	
	  ;% helicopter_dag_2_P.HILInitialize_DOPEnter
	  section.data(20).logicalSrcIdx = 102;
	  section.data(20).dtTransOffset = 19;
	
	  ;% helicopter_dag_2_P.HILInitialize_DOStart
	  section.data(21).logicalSrcIdx = 103;
	  section.data(21).dtTransOffset = 20;
	
	  ;% helicopter_dag_2_P.HILInitialize_DOEnter
	  section.data(22).logicalSrcIdx = 104;
	  section.data(22).dtTransOffset = 21;
	
	  ;% helicopter_dag_2_P.HILInitialize_DOReset
	  section.data(23).logicalSrcIdx = 105;
	  section.data(23).dtTransOffset = 22;
	
	  ;% helicopter_dag_2_P.HILInitialize_EIPStart
	  section.data(24).logicalSrcIdx = 106;
	  section.data(24).dtTransOffset = 23;
	
	  ;% helicopter_dag_2_P.HILInitialize_EIPEnter
	  section.data(25).logicalSrcIdx = 107;
	  section.data(25).dtTransOffset = 24;
	
	  ;% helicopter_dag_2_P.HILInitialize_EIStart
	  section.data(26).logicalSrcIdx = 108;
	  section.data(26).dtTransOffset = 25;
	
	  ;% helicopter_dag_2_P.HILInitialize_EIEnter
	  section.data(27).logicalSrcIdx = 109;
	  section.data(27).dtTransOffset = 26;
	
	  ;% helicopter_dag_2_P.HILInitialize_POPStart
	  section.data(28).logicalSrcIdx = 110;
	  section.data(28).dtTransOffset = 27;
	
	  ;% helicopter_dag_2_P.HILInitialize_POPEnter
	  section.data(29).logicalSrcIdx = 111;
	  section.data(29).dtTransOffset = 28;
	
	  ;% helicopter_dag_2_P.HILInitialize_POStart
	  section.data(30).logicalSrcIdx = 112;
	  section.data(30).dtTransOffset = 29;
	
	  ;% helicopter_dag_2_P.HILInitialize_POEnter
	  section.data(31).logicalSrcIdx = 113;
	  section.data(31).dtTransOffset = 30;
	
	  ;% helicopter_dag_2_P.HILInitialize_POReset
	  section.data(32).logicalSrcIdx = 114;
	  section.data(32).dtTransOffset = 31;
	
	  ;% helicopter_dag_2_P.HILInitialize_OOReset
	  section.data(33).logicalSrcIdx = 115;
	  section.data(33).dtTransOffset = 32;
	
	  ;% helicopter_dag_2_P.HILInitialize_DOFinal
	  section.data(34).logicalSrcIdx = 116;
	  section.data(34).dtTransOffset = 33;
	
	  ;% helicopter_dag_2_P.HILInitialize_DOInitial
	  section.data(35).logicalSrcIdx = 117;
	  section.data(35).dtTransOffset = 34;
	
	  ;% helicopter_dag_2_P.HILReadEncoderTimebase_Active
	  section.data(36).logicalSrcIdx = 118;
	  section.data(36).dtTransOffset = 35;
	
	  ;% helicopter_dag_2_P.HILWriteAnalog_Active
	  section.data(37).logicalSrcIdx = 119;
	  section.data(37).dtTransOffset = 36;
	
	  ;% helicopter_dag_2_P.GameController_AutoCenter
	  section.data(38).logicalSrcIdx = 120;
	  section.data(38).dtTransOffset = 37;
	
	  ;% helicopter_dag_2_P.GameController_Enabled
	  section.data(39).logicalSrcIdx = 121;
	  section.data(39).dtTransOffset = 38;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(7) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_P.HILReadEncoderTimebase_Overflow
	  section.data(1).logicalSrcIdx = 122;
	  section.data(1).dtTransOffset = 0;
	
	  ;% helicopter_dag_2_P.GameController_ControllerNumber
	  section.data(2).logicalSrcIdx = 123;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(8) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (parameter)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    paramMap.nTotData = nTotData;
    


  ;%**************************
  ;% Create Block Output Map *
  ;%**************************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 1;
    sectIdxOffset = 0;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc sigMap
    ;%
    sigMap.nSections           = nTotSects;
    sigMap.sectIdxOffset       = sectIdxOffset;
      sigMap.sections(nTotSects) = dumSection; %prealloc
    sigMap.nTotData            = -1;
    
    ;%
    ;% Auto data (helicopter_dag_2_B)
    ;%
      section.nData     = 21;
      section.data(21)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_B.FromWorkspace
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% helicopter_dag_2_B.TravelCounttorad
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% helicopter_dag_2_B.Gain
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 2;
	
	  ;% helicopter_dag_2_B.Gain_d
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 3;
	
	  ;% helicopter_dag_2_B.PitchCounttorad
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 4;
	
	  ;% helicopter_dag_2_B.Gain_i
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 5;
	
	  ;% helicopter_dag_2_B.Gain_b
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 6;
	
	  ;% helicopter_dag_2_B.ElevationCounttorad
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 7;
	
	  ;% helicopter_dag_2_B.Gain_e
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 8;
	
	  ;% helicopter_dag_2_B.Sum
	  section.data(10).logicalSrcIdx = 9;
	  section.data(10).dtTransOffset = 9;
	
	  ;% helicopter_dag_2_B.Gain_dg
	  section.data(11).logicalSrcIdx = 10;
	  section.data(11).dtTransOffset = 10;
	
	  ;% helicopter_dag_2_B.Gain1
	  section.data(12).logicalSrcIdx = 11;
	  section.data(12).dtTransOffset = 11;
	
	  ;% helicopter_dag_2_B.TmpSignalConversionAtToFileInpo
	  section.data(13).logicalSrcIdx = 12;
	  section.data(13).dtTransOffset = 17;
	
	  ;% helicopter_dag_2_B.RateTransitiony
	  section.data(14).logicalSrcIdx = 13;
	  section.data(14).dtTransOffset = 19;
	
	  ;% helicopter_dag_2_B.Joystick_gain_y
	  section.data(15).logicalSrcIdx = 14;
	  section.data(15).dtTransOffset = 20;
	
	  ;% helicopter_dag_2_B.Gain_l
	  section.data(16).logicalSrcIdx = 15;
	  section.data(16).dtTransOffset = 21;
	
	  ;% helicopter_dag_2_B.BackmotorSaturation
	  section.data(17).logicalSrcIdx = 16;
	  section.data(17).dtTransOffset = 22;
	
	  ;% helicopter_dag_2_B.FrontmotorSaturation
	  section.data(18).logicalSrcIdx = 17;
	  section.data(18).dtTransOffset = 23;
	
	  ;% helicopter_dag_2_B.RateTransitionx
	  section.data(19).logicalSrcIdx = 18;
	  section.data(19).dtTransOffset = 24;
	
	  ;% helicopter_dag_2_B.Joystick_gain_x
	  section.data(20).logicalSrcIdx = 19;
	  section.data(20).dtTransOffset = 25;
	
	  ;% helicopter_dag_2_B.In1
	  section.data(21).logicalSrcIdx = 20;
	  section.data(21).dtTransOffset = 26;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(1) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (signal)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    sigMap.nTotData = nTotData;
    


  ;%*******************
  ;% Create DWork Map *
  ;%*******************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 9;
    sectIdxOffset = 1;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc dworkMap
    ;%
    dworkMap.nSections           = nTotSects;
    dworkMap.sectIdxOffset       = sectIdxOffset;
      dworkMap.sections(nTotSects) = dumSection; %prealloc
    dworkMap.nTotData            = -1;
    
    ;%
    ;% Auto data (helicopter_dag_2_DW)
    ;%
      section.nData     = 15;
      section.data(15)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_DW.HILInitialize_AIMinimums
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% helicopter_dag_2_DW.HILInitialize_AIMaximums
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 8;
	
	  ;% helicopter_dag_2_DW.HILInitialize_AOMinimums
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 16;
	
	  ;% helicopter_dag_2_DW.HILInitialize_AOMaximums
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 24;
	
	  ;% helicopter_dag_2_DW.HILInitialize_AOVoltages
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 32;
	
	  ;% helicopter_dag_2_DW.HILInitialize_FilterFrequency
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 40;
	
	  ;% helicopter_dag_2_DW.HILInitialize_POSortedFreqs
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 48;
	
	  ;% helicopter_dag_2_DW.HILInitialize_POValues
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 56;
	
	  ;% helicopter_dag_2_DW.RateTransitiony_Buffer0
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 64;
	
	  ;% helicopter_dag_2_DW.TimeStampA
	  section.data(10).logicalSrcIdx = 9;
	  section.data(10).dtTransOffset = 65;
	
	  ;% helicopter_dag_2_DW.LastUAtTimeA
	  section.data(11).logicalSrcIdx = 10;
	  section.data(11).dtTransOffset = 66;
	
	  ;% helicopter_dag_2_DW.TimeStampB
	  section.data(12).logicalSrcIdx = 11;
	  section.data(12).dtTransOffset = 67;
	
	  ;% helicopter_dag_2_DW.LastUAtTimeB
	  section.data(13).logicalSrcIdx = 12;
	  section.data(13).dtTransOffset = 68;
	
	  ;% helicopter_dag_2_DW.HILWriteAnalog_Buffer
	  section.data(14).logicalSrcIdx = 13;
	  section.data(14).dtTransOffset = 69;
	
	  ;% helicopter_dag_2_DW.RateTransitionx_Buffer0
	  section.data(15).logicalSrcIdx = 14;
	  section.data(15).dtTransOffset = 71;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(1) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_DW.GameController_Controller
	  section.data(1).logicalSrcIdx = 15;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(2) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_DW.HILInitialize_Card
	  section.data(1).logicalSrcIdx = 16;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(3) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_DW.HILReadEncoderTimebase_Task
	  section.data(1).logicalSrcIdx = 17;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(4) = section;
      clear section
      
      section.nData     = 15;
      section.data(15)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_DW.FromWorkspace_PWORK.TimePtr
	  section.data(1).logicalSrcIdx = 18;
	  section.data(1).dtTransOffset = 0;
	
	  ;% helicopter_dag_2_DW.ToFile_PWORK.FilePtr
	  section.data(2).logicalSrcIdx = 19;
	  section.data(2).dtTransOffset = 1;
	
	  ;% helicopter_dag_2_DW.lambdavsp_c_PWORK.LoggedData
	  section.data(3).logicalSrcIdx = 20;
	  section.data(3).dtTransOffset = 2;
	
	  ;% helicopter_dag_2_DW.ElevationScopedegs_PWORK.LoggedData
	  section.data(4).logicalSrcIdx = 21;
	  section.data(4).dtTransOffset = 3;
	
	  ;% helicopter_dag_2_DW.ElevationScopedeg_PWORK.LoggedData
	  section.data(5).logicalSrcIdx = 22;
	  section.data(5).dtTransOffset = 4;
	
	  ;% helicopter_dag_2_DW.PitchScopedeg_PWORK.LoggedData
	  section.data(6).logicalSrcIdx = 23;
	  section.data(6).dtTransOffset = 5;
	
	  ;% helicopter_dag_2_DW.PtichrateScopedegs_PWORK.LoggedData
	  section.data(7).logicalSrcIdx = 24;
	  section.data(7).dtTransOffset = 6;
	
	  ;% helicopter_dag_2_DW.PtichrateScopedegs1_PWORK.LoggedData
	  section.data(8).logicalSrcIdx = 25;
	  section.data(8).dtTransOffset = 7;
	
	  ;% helicopter_dag_2_DW.TravelrateScopedegs_PWORK.LoggedData
	  section.data(9).logicalSrcIdx = 26;
	  section.data(9).dtTransOffset = 8;
	
	  ;% helicopter_dag_2_DW.TravelScopedeg_PWORK.LoggedData
	  section.data(10).logicalSrcIdx = 27;
	  section.data(10).dtTransOffset = 9;
	
	  ;% helicopter_dag_2_DW.Backmotor_PWORK.LoggedData
	  section.data(11).logicalSrcIdx = 28;
	  section.data(11).dtTransOffset = 10;
	
	  ;% helicopter_dag_2_DW.Frontmotor_PWORK.LoggedData
	  section.data(12).logicalSrcIdx = 29;
	  section.data(12).dtTransOffset = 11;
	
	  ;% helicopter_dag_2_DW.HILWriteAnalog_PWORK
	  section.data(13).logicalSrcIdx = 30;
	  section.data(13).dtTransOffset = 12;
	
	  ;% helicopter_dag_2_DW.XScope_PWORK.LoggedData
	  section.data(14).logicalSrcIdx = 31;
	  section.data(14).dtTransOffset = 13;
	
	  ;% helicopter_dag_2_DW.YScope_PWORK.LoggedData
	  section.data(15).logicalSrcIdx = 32;
	  section.data(15).dtTransOffset = 14;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(5) = section;
      clear section
      
      section.nData     = 7;
      section.data(7)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_DW.HILInitialize_ClockModes
	  section.data(1).logicalSrcIdx = 33;
	  section.data(1).dtTransOffset = 0;
	
	  ;% helicopter_dag_2_DW.HILInitialize_QuadratureModes
	  section.data(2).logicalSrcIdx = 34;
	  section.data(2).dtTransOffset = 3;
	
	  ;% helicopter_dag_2_DW.HILInitialize_InitialEICounts
	  section.data(3).logicalSrcIdx = 35;
	  section.data(3).dtTransOffset = 11;
	
	  ;% helicopter_dag_2_DW.HILInitialize_POModeValues
	  section.data(4).logicalSrcIdx = 36;
	  section.data(4).dtTransOffset = 19;
	
	  ;% helicopter_dag_2_DW.HILInitialize_POAlignValues
	  section.data(5).logicalSrcIdx = 37;
	  section.data(5).dtTransOffset = 27;
	
	  ;% helicopter_dag_2_DW.HILInitialize_POPolarityVals
	  section.data(6).logicalSrcIdx = 38;
	  section.data(6).dtTransOffset = 35;
	
	  ;% helicopter_dag_2_DW.HILReadEncoderTimebase_Buffer
	  section.data(7).logicalSrcIdx = 39;
	  section.data(7).dtTransOffset = 43;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(6) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_DW.HILInitialize_POSortedChans
	  section.data(1).logicalSrcIdx = 40;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(7) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_DW.FromWorkspace_IWORK.PrevIndex
	  section.data(1).logicalSrcIdx = 41;
	  section.data(1).dtTransOffset = 0;
	
	  ;% helicopter_dag_2_DW.ToFile_IWORK.Count
	  section.data(2).logicalSrcIdx = 42;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(8) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% helicopter_dag_2_DW.If_ActiveSubsystem
	  section.data(1).logicalSrcIdx = 43;
	  section.data(1).dtTransOffset = 0;
	
	  ;% helicopter_dag_2_DW.IfActionSubsystem_SubsysRanBC
	  section.data(2).logicalSrcIdx = 44;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(9) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (dwork)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    dworkMap.nTotData = nTotData;
    


  ;%
  ;% Add individual maps to base struct.
  ;%

  targMap.paramMap  = paramMap;    
  targMap.signalMap = sigMap;
  targMap.dworkMap  = dworkMap;
  
  ;%
  ;% Add checksums to base struct.
  ;%


  targMap.checksum0 = 266776760;
  targMap.checksum1 = 680672581;
  targMap.checksum2 = 2076109412;
  targMap.checksum3 = 2763510147;

