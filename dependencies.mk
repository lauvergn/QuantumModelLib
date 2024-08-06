#===============================================
adiachannels_basis_m = $(OBJ_DIR)/Basis_m.o
adiachannels_makehinact_m = $(OBJ_DIR)/MakeHinact_m.o
model_m = $(OBJ_DIR)/Model_m.o
irc_m = $(OBJ_DIR)/IRC_m.o
opt_m = $(OBJ_DIR)/Opt_m.o
qml_bottleneck_m = $(OBJ_DIR)/Bottleneck_m.o
qml_buck_m = $(OBJ_DIR)/Buck_m.o
qml_ch5_m = $(OBJ_DIR)/CH5_m.o
qml_cnh_murrell_m = $(OBJ_DIR)/CNH_Murrell_m.o
qml_clh2p_botschwina_m = $(OBJ_DIR)/ClH2p_Botschwina_m.o
qml_clh2p_m = $(OBJ_DIR)/ClH2p_m.o
qml_empty_m = $(OBJ_DIR)/Empty_m.o
qml_h2nsi_m = $(OBJ_DIR)/H2NSi_m.o
qml_h2o_m = $(OBJ_DIR)/H2O_m.o
qml_h2sin_m = $(OBJ_DIR)/H2SiN_m.o
qml_h2_m = $(OBJ_DIR)/H2_m.o
qml_h3_m = $(OBJ_DIR)/H3_m.o
qml_hnnhp_m = $(OBJ_DIR)/HNNHp_m.o
qml_hno3_m = $(OBJ_DIR)/HNO3_m.o
qml_hono_m = $(OBJ_DIR)/HONO_m.o
qml_hoo_dmbe_m = $(OBJ_DIR)/HOO_DMBE_m.o
qml_henonheiles_m = $(OBJ_DIR)/HenonHeiles_m.o
qml_linearhbond_m = $(OBJ_DIR)/LinearHBond_m.o
qml_morse_m = $(OBJ_DIR)/Morse_m.o
qml_no3_m = $(OBJ_DIR)/NO3_m.o
qml_onedsoc_1s1t_m = $(OBJ_DIR)/OneDSOC_1S1T_m.o
qml_onedsoc_2s1t_m = $(OBJ_DIR)/OneDSOC_2S1T_m.o
qml_oned_photons_m = $(OBJ_DIR)/OneD_Photons_m.o
qml_ph4jo_m = $(OBJ_DIR)/PH4Jo_m.o
qml_ph4_m = $(OBJ_DIR)/PH4_m.o
qml_psb3_m = $(OBJ_DIR)/PSB3_m.o
qml_phenol_m = $(OBJ_DIR)/Phenol_m.o
qml_poly1d_m = $(OBJ_DIR)/Poly1D_m.o
qml_retinal_jpcb2000_m = $(OBJ_DIR)/Retinal_JPCB2000_m.o
qml_sigmoid_m = $(OBJ_DIR)/Sigmoid_m.o
qml_template_m = $(OBJ_DIR)/Template_m.o
qml_test_m = $(OBJ_DIR)/Test_m.o
qml_tully_m = $(OBJ_DIR)/Tully_m.o
qml_twod_mullerbrown_m = $(OBJ_DIR)/TwoD_MullerBrown_m.o
qml_twod_rjdi2014_m = $(OBJ_DIR)/TwoD_RJDI2014_m.o
qml_twod_valahu2022_m = $(OBJ_DIR)/TwoD_Valahu2022_m.o
qml_twod_m = $(OBJ_DIR)/TwoD_m.o
qml_uracil_m = $(OBJ_DIR)/Uracil_m.o
qml_vibronic_m = $(OBJ_DIR)/Vibronic_m.o
qmllib_finitediff_m = $(OBJ_DIR)/FiniteDiff_m.o
qmllib_utillib_m = $(OBJ_DIR)/UtilLib_m.o
#===============================================
$(OBJ_DIR)/Basis_m.o : \
          $(qdutil_numparameters_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/MakeHinact_m.o : \
          $(qdutil_numparameters_m) \
          $(qdutil_m) \
          $(addnsvm_m) \
          $(model_m) \
          $(adiachannels_basis_m)
$(OBJ_DIR)/Model_driver.o : \
          $(qdutil_numparameters_m) \
          $(model_m) \
          $(opt_m) \
          $(addnsvm_m) \
          $(qdutil_m)
$(OBJ_DIR)/Model_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(adiachannels_basis_m) \
          $(qdutil_m) \
          $(qmllib_utillib_m) \
          $(qml_template_m) \
          $(qml_test_m) \
          $(qml_morse_m) \
          $(qml_poly1d_m) \
          $(qml_h2_m) \
          $(qml_henonheiles_m) \
          $(qml_tully_m) \
          $(qml_psb3_m) \
          $(qml_retinal_jpcb2000_m) \
          $(qml_hono_m) \
          $(qml_hnnhp_m) \
          $(qml_h2sin_m) \
          $(qml_h2nsi_m) \
          $(qml_h2o_m) \
          $(qml_clh2p_m) \
          $(qml_clh2p_botschwina_m) \
          $(qml_bottleneck_m) \
          $(qml_hno3_m) \
          $(qml_no3_m) \
          $(qml_ch5_m) \
          $(qml_ph4jo_m) \
          $(qml_ph4_m) \
          $(qml_hoo_dmbe_m) \
          $(qml_h3_m) \
          $(qml_cnh_murrell_m) \
          $(qml_onedsoc_1s1t_m) \
          $(qml_onedsoc_2s1t_m) \
          $(qml_linearhbond_m) \
          $(qml_twod_mullerbrown_m) \
          $(qml_buck_m) \
          $(qml_phenol_m) \
          $(qml_sigmoid_m) \
          $(qml_twod_m) \
          $(qml_twod_rjdi2014_m) \
          $(qml_twod_valahu2022_m) \
          $(qml_vibronic_m) \
          $(qml_uracil_m) \
          $(qml_oned_photons_m) \
          $(addnsvm_m) \
          $(qmllib_finitediff_m) \
          $(qdutil_test_m)
$(OBJ_DIR)/IRC_m.o : \
          $(qdutil_numparameters_m) \
          $(opt_m) \
          $(qmllib_utillib_m) \
          $(qdutil_m) \
          $(model_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Opt_m.o : \
          $(qdutil_numparameters_m) \
          $(qdutil_m) \
          $(model_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Bottleneck_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Buck_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(addnsvm_m)
$(OBJ_DIR)/CH5_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(qmllib_utillib_m) \
          $(addnsvm_m)
$(OBJ_DIR)/CNH_Murrell_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/ClH2p_Botschwina_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(addnsvm_m)
$(OBJ_DIR)/ClH2p_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(qmllib_utillib_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Empty_m.o : \
          $(qdutil_numparameters_m) \
          $(addnsvm_m) \
          $(qdutil_m) \
          $(qmllib_utillib_m)
$(OBJ_DIR)/H2NSi_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(qmllib_utillib_m) \
          $(addnsvm_m)
$(OBJ_DIR)/H2O_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(addnsvm_m)
$(OBJ_DIR)/H2SiN_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(qmllib_utillib_m) \
          $(addnsvm_m)
$(OBJ_DIR)/H2_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(addnsvm_m)
$(OBJ_DIR)/H3_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/HNNHp_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(qmllib_utillib_m) \
          $(addnsvm_m)
$(OBJ_DIR)/HNO3_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(qmllib_utillib_m) \
          $(addnsvm_m)
$(OBJ_DIR)/HONO_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/HOO_DMBE_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/HenonHeiles_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qml_morse_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/LinearHBond_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qml_morse_m) \
          $(qml_buck_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Morse_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(addnsvm_m)
$(OBJ_DIR)/NO3_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/OneDSOC_1S1T_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/OneDSOC_2S1T_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/OneD_Photons_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/PH4Jo_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(qmllib_utillib_m) \
          $(addnsvm_m)
$(OBJ_DIR)/PH4_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(qmllib_utillib_m) \
          $(addnsvm_m)
$(OBJ_DIR)/PSB3_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Phenol_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qml_morse_m) \
          $(qml_sigmoid_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Poly1D_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Retinal_JPCB2000_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Sigmoid_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Template_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qml_morse_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Test_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Tully_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/TwoD_MullerBrown_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/TwoD_RJDI2014_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/TwoD_Valahu2022_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/TwoD_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Uracil_m.o : \
          $(qdutil_numparameters_m) \
          $(qml_empty_m) \
          $(qdutil_m) \
          $(qmllib_utillib_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Vibronic_m.o : \
          $(qdutil_numparameters_m) \
          $(addnsvm_m) \
          $(qml_empty_m) \
          $(qdutil_m)
$(OBJ_DIR)/FiniteDiff_m.o : \
          $(qdutil_numparameters_m) \
          $(addnsvm_m)
$(OBJ_DIR)/UtilLib_m.o : \
          $(qdutil_numparameters_m) \
          $(qdutil_m)
