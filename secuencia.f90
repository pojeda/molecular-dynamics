SUBROUTINE SECUENCIAS

! EPSILONS AND SIGMAS USED FOR THE LENNARD JONES POTENTIAL
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS

	IF(SECU.EQ.1) THEN
		DO I=1,NUM_RES
			DO J=1,NUM_RES
				SIGMA_CONST(I,J)=6.5
			ENDDO
		ENDDO
		
		DO I=1,4
			DO J=1,4
				EPS_CONST(I,J)=40.
			ENDDO
		ENDDO
	ENDIF

	IF(SECU.EQ.2.OR.SECU.EQ.3.OR.SECU.EQ.4) THEN
	DO I=1,NUM_RES
		DO J=1,NUM_RES
			SIGMA_CONST(I,J)=6.5
		ENDDO
	ENDDO

	eps_const(1,1)=40.0
	eps_const(1,2)=30.0
	eps_const(1,3)=20.0
	eps_const(1,4)=17.0
	eps_const(2,1)=30.0
	eps_const(2,2)=25.0
	eps_const(2,3)=13.0
	eps_const(2,4)=10.0
	eps_const(3,1)=20.0
	eps_const(3,2)=13.0
	eps_const(3,3)=5.0
	eps_const(3,4)=2.0
	eps_const(4,1)=17.0
	eps_const(4,2)=10.0
	eps_const(4,3)=2.0
	eps_const(4,4)=1.0
	ENDIF

	IF(SECU.NE.1.AND.SECU.NE.2.AND.SECU.NE.3.AND.SECU.NE.4) THEN
	sigma_const( 1, 3)=  6.459669370163
	sigma_const( 1, 4)=  8.968331577040
	sigma_const( 1, 5)= 10.723693904198
	sigma_const( 1, 6)= 12.150344895811
	sigma_const( 1, 7)= 13.566126020422
	sigma_const( 1, 8)= 14.351283961647
	sigma_const( 1, 9)= 14.550255764072
	sigma_const( 1,10)= 15.151319369087
	sigma_const( 1,11)= 15.290811357634
	sigma_const( 1,12)= 13.699298422737
	sigma_const( 1,13)= 12.636084696811
	sigma_const( 1,14)= 12.400390700090
	sigma_const( 1,15)= 12.443398211504
	sigma_const( 1,16)= 11.637801304320
	sigma_const( 1,17)= 10.457029380276
	sigma_const( 1,18)=  8.993109775970
	sigma_const( 1,19)=  7.164499962542
	sigma_const( 1,20)=  6.269935597371
	sigma_const( 1,21)=  6.397981932181
	sigma_const( 1,22)=  7.851560685608
	sigma_const( 1,23)=  9.261807574385
	sigma_const( 1,24)=  8.811123738126
	sigma_const( 1,25)=  9.199218744045
	sigma_const( 1,26)=  9.381221708965
	sigma_const( 1,27)=  9.433085597544
	sigma_const( 1,28)=  8.347798886489
	sigma_const( 1,29)=  7.062508164594
	sigma_const( 1,30)=  5.797551109593
	sigma_const( 2, 4)=  6.319624810247
	sigma_const( 2, 5)=  8.802203228910
	sigma_const( 2, 6)= 10.844334511405
	sigma_const( 2, 7)= 12.260227568205
	sigma_const( 2, 8)= 12.349604216665
	sigma_const( 2, 9)= 11.932293693109
	sigma_const( 2,10)= 12.218542736091
	sigma_const( 2,11)= 12.466527913970
	sigma_const( 2,12)= 10.943458950271
	sigma_const( 2,13)=  9.526843173066
	sigma_const( 2,14)=  9.041794182378
	sigma_const( 2,15)=  9.184422564600
	sigma_const( 2,16)=  8.811659462832
	sigma_const( 2,17)=  8.317758558443
	sigma_const( 2,18)=  7.716511804375
	sigma_const( 2,19)=  6.824755852443
	sigma_const( 2,20)=  6.510762004872
	sigma_const( 2,21)=  6.454990333068
	sigma_const( 2,22)=  6.876631192286
	sigma_const( 2,23)=  7.280410660536
	sigma_const( 2,24)=  6.533474831376
	sigma_const( 2,25)=  6.531401251613
	sigma_const( 2,26)=  6.577633693262
	sigma_const( 2,27)=  6.760464468418
	sigma_const( 2,28)=  6.428533427053
	sigma_const( 2,29)=  6.310835500392
	sigma_const( 2,30)=  6.575462291344
	sigma_const( 3, 5)=  6.380866259812
	sigma_const( 3, 6)=  8.797946780773
	sigma_const( 3, 7)= 10.151897103370
	sigma_const( 3, 8)=  9.659760458176
	sigma_const( 3, 9)=  8.776824381591
	sigma_const( 3,10)=  8.901282347310
	sigma_const( 3,11)=  9.374365752930
	sigma_const( 3,12)=  8.320974917384
	sigma_const( 3,13)=  7.023201358528
	sigma_const( 3,14)=  6.340084535521
	sigma_const( 3,15)=  6.366767853318
	sigma_const( 3,16)=  6.333349088626
	sigma_const( 3,17)=  6.313414139007
	sigma_const( 3,18)=  6.412645146877
	sigma_const( 3,19)=  6.442411703145
	sigma_const( 3,20)=  6.719391670973
	sigma_const( 3,21)=  6.859455926827
	sigma_const( 3,22)=  6.697644924825
	sigma_const( 3,23)=  6.437983960761
	sigma_const( 3,24)=  6.315914152138
	sigma_const( 3,25)=  6.381750823318
	sigma_const( 3,26)=  6.365814334970
	sigma_const( 3,27)=  6.328136740098
	sigma_const( 3,28)=  6.555942774702
	sigma_const( 3,29)=  7.107034336886
	sigma_const( 3,30)=  8.098306477303
	sigma_const( 4, 6)=  6.151175538571
	sigma_const( 4, 7)=  7.873493456252
	sigma_const( 4, 8)=  7.321239610154
	sigma_const( 4, 9)=  6.345136818351
	sigma_const( 4,10)=  6.304611847174
	sigma_const( 4,11)=  6.611680123638
	sigma_const( 4,12)=  6.259524772881
	sigma_const( 4,13)=  6.224454932848
	sigma_const( 4,14)=  6.264247073130
	sigma_const( 4,15)=  6.405203969803
	sigma_const( 4,16)=  6.616126184956
	sigma_const( 4,17)=  6.420834868307
	sigma_const( 4,18)=  6.333194867333
	sigma_const( 4,19)=  6.463995711388
	sigma_const( 4,20)=  6.652969749751
	sigma_const( 4,21)=  6.811878170537
	sigma_const( 4,22)=  6.458037277209
	sigma_const( 4,23)=  6.348492209098
	sigma_const( 4,24)=  7.604255348886
	sigma_const( 4,25)=  8.490583953948
	sigma_const( 4,26)=  8.953678963286
	sigma_const( 4,27)=  8.945144490386
	sigma_const( 4,28)=  9.205800371502
	sigma_const( 4,29)=  9.520446274085
	sigma_const( 4,30)= 10.206060280676
	sigma_const( 5, 7)=  6.078751750417
	sigma_const( 5, 8)=  6.497299617475
	sigma_const( 5, 9)=  6.437747223287
	sigma_const( 5,10)=  6.628829194819
	sigma_const( 5,11)=  6.389904675702
	sigma_const( 5,12)=  6.706277234792
	sigma_const( 5,13)=  8.178849100939
	sigma_const( 5,14)=  8.978248006919
	sigma_const( 5,15)=  9.155509418870
	sigma_const( 5,16)=  9.107766798760
	sigma_const( 5,17)=  8.207472767351
	sigma_const( 5,18)=  7.210332482415
	sigma_const( 5,19)=  6.670938686117
	sigma_const( 5,20)=  6.341060612846
	sigma_const( 5,21)=  6.579771299576
	sigma_const( 5,22)=  6.671086307104
	sigma_const( 5,23)=  7.507421059231
	sigma_const( 5,24)=  9.818347166576
	sigma_const( 5,25)= 11.412023874332
	sigma_const( 5,26)= 12.242165612577
	sigma_const( 5,27)= 12.245273507564
	sigma_const( 5,28)= 12.169555368793
	sigma_const( 5,29)= 11.884925976109
	sigma_const( 5,30)= 11.872983141730
	sigma_const( 6, 8)=  5.366622511954
	sigma_const( 6, 9)=  6.984267341906
	sigma_const( 6,10)=  8.361422843704
	sigma_const( 6,11)=  8.773503585005
	sigma_const( 6,12)=  9.821413929552
	sigma_const( 6,13)= 11.445709855497
	sigma_const( 6,14)= 11.878332059913
	sigma_const( 6,15)= 11.358944688075
	sigma_const( 6,16)= 10.471526629010
	sigma_const( 6,17)=  8.667614131056
	sigma_const( 6,18)=  6.904783848307
	sigma_const( 6,19)=  6.293697421927
	sigma_const( 6,20)=  6.577981618759
	sigma_const( 6,21)=  7.957191722638
	sigma_const( 6,22)=  9.140134567523
	sigma_const( 6,23)= 10.656575166754
	sigma_const( 6,24)= 13.002622149047
	sigma_const( 6,25)= 14.491444135288
	sigma_const( 6,26)= 14.970465248451
	sigma_const( 6,27)= 14.423880582239
	sigma_const( 6,28)= 13.663596271154
	sigma_const( 6,29)= 12.694737370327
	sigma_const( 6,30)= 12.168262776453
	sigma_const( 7, 9)=  6.382910709869
	sigma_const( 7,10)=  8.931762495764
	sigma_const( 7,11)= 10.432688007702
	sigma_const( 7,12)= 12.122815471889
	sigma_const( 7,13)= 13.324191660825
	sigma_const( 7,14)= 12.976394321031
	sigma_const( 7,15)= 11.584180955448
	sigma_const( 7,16)=  9.928237735210
	sigma_const( 7,17)=  7.656545071927
	sigma_const( 7,18)=  6.096908526878
	sigma_const( 7,19)=  6.733100055952
	sigma_const( 7,20)=  8.400688728949
	sigma_const( 7,21)= 10.603926095928
	sigma_const( 7,22)= 12.093920674729
	sigma_const( 7,23)= 13.416815429001
	sigma_const( 7,24)= 15.323398551909
	sigma_const( 7,25)= 16.236055146175
	sigma_const( 7,26)= 16.048495065853
	sigma_const( 7,27)= 14.807409112941
	sigma_const( 7,28)= 13.571063526471
	sigma_const( 7,29)= 12.394044588468
	sigma_const( 7,30)= 12.062193718164
	sigma_const( 8,10)=  6.494267108035
	sigma_const( 8,11)=  8.857710982033
	sigma_const( 8,12)= 11.061486901001
	sigma_const( 8,13)= 11.856275739425
	sigma_const( 8,14)= 10.974780976643
	sigma_const( 8,15)=  9.193564909731
	sigma_const( 8,16)=  7.684603683205
	sigma_const( 8,17)=  6.240094659245
	sigma_const( 8,18)=  6.326182668396
	sigma_const( 8,19)=  8.211168855331
	sigma_const( 8,20)= 10.215954976249
	sigma_const( 8,21)= 12.087676336342
	sigma_const( 8,22)= 12.902725153348
	sigma_const( 8,23)= 13.397616947782
	sigma_const( 8,24)= 14.874785712637
	sigma_const( 8,25)= 15.282304641665
	sigma_const( 8,26)= 14.748264759272
	sigma_const( 8,27)= 13.342462878462
	sigma_const( 8,28)= 12.512991062607
	sigma_const( 8,29)= 12.094125343237
	sigma_const( 8,30)= 12.689425501727
	sigma_const( 9,11)=  6.360977652835
	sigma_const( 9,12)=  8.793871605818
	sigma_const( 9,13)=  9.266242294784
	sigma_const( 9,14)=  8.216572055287
	sigma_const( 9,15)=  6.647664371572
	sigma_const( 9,16)=  6.219123846093
	sigma_const( 9,17)=  6.388738413441
	sigma_const( 9,18)=  7.702187291514
	sigma_const( 9,19)=  9.712705676576
	sigma_const( 9,20)= 11.332063305837
	sigma_const( 9,21)= 12.488024674914
	sigma_const( 9,22)= 12.448003077219
	sigma_const( 9,23)= 12.081850853235
	sigma_const( 9,24)= 13.240867470611
	sigma_const( 9,25)= 13.396452042838
	sigma_const( 9,26)= 12.914420974393
	sigma_const( 9,27)= 11.853366338188
	sigma_const( 9,28)= 11.862105967462
	sigma_const( 9,29)= 12.360217866052
	sigma_const( 9,30)= 13.634852169628
	sigma_const(10,12)=  6.106565883180
	sigma_const(10,13)=  6.807795257211
	sigma_const(10,14)=  6.459570573267
	sigma_const(10,15)=  6.227829940749
	sigma_const(10,16)=  7.470796225420
	sigma_const(10,17)=  8.655910733592
	sigma_const(10,18)= 10.057534398158
	sigma_const(10,19)= 11.528604662852
	sigma_const(10,20)= 12.384919561569
	sigma_const(10,21)= 12.651308109526
	sigma_const(10,22)= 11.711434205072
	sigma_const(10,23)= 10.592333524565
	sigma_const(10,24)= 11.766335255619
	sigma_const(10,25)= 12.122644734767
	sigma_const(10,26)= 12.215154336206
	sigma_const(10,27)= 11.978673907803
	sigma_const(10,28)= 12.871807791052
	sigma_const(10,29)= 13.955374342591
	sigma_const(10,30)= 15.394375905921
	sigma_const(11,13)=  5.506008657634
	sigma_const(11,14)=  6.889047991457
	sigma_const(11,15)=  8.145400416160
	sigma_const(11,16)=  9.970023237084
	sigma_const(11,17)= 11.049804715564
	sigma_const(11,18)= 11.862994393043
	sigma_const(11,19)= 12.497988428686
	sigma_const(11,20)= 12.447561372152
	sigma_const(11,21)= 11.858713715914
	sigma_const(11,22)= 10.228347913848
	sigma_const(11,23)=  8.813411985282
	sigma_const(11,24)= 10.512604253117
	sigma_const(11,25)= 11.630190300263
	sigma_const(11,26)= 12.639352479369
	sigma_const(11,27)= 13.252382591346
	sigma_const(11,28)= 14.521053522633
	sigma_const(11,29)= 15.597742276955
	sigma_const(11,30)= 16.707165077402
	sigma_const(12,14)=  6.191971366680
	sigma_const(12,15)=  8.551014523799
	sigma_const(12,16)= 10.747515217616
	sigma_const(12,17)= 11.905911552677
	sigma_const(12,18)= 12.456321307433
	sigma_const(12,19)= 12.484823935672
	sigma_const(12,20)= 11.748149736055
	sigma_const(12,21)= 10.354512687159
	sigma_const(12,22)=  7.986055837725
	sigma_const(12,23)=  5.900408405208
	sigma_const(12,24)=  7.621763050335
	sigma_const(12,25)=  9.247851725094
	sigma_const(12,26)= 11.007390002276
	sigma_const(12,27)= 12.394307903713
	sigma_const(12,28)= 14.028071172711
	sigma_const(12,29)= 15.232594278591
	sigma_const(12,30)= 16.213942879658
	sigma_const(13,15)=  6.454495379195
	sigma_const(13,16)=  9.174073795385
	sigma_const(13,17)= 10.994810702335
	sigma_const(13,18)= 12.139561778575
	sigma_const(13,19)= 12.554775854292
	sigma_const(13,20)= 12.188867095672
	sigma_const(13,21)= 10.951972290842
	sigma_const(13,22)=  8.602602116000
	sigma_const(13,23)=  5.848548907096
	sigma_const(13,24)=  5.914136011610
	sigma_const(13,25)=  6.472512752319
	sigma_const(13,26)=  7.972107745299
	sigma_const(13,27)=  9.627730681977
	sigma_const(13,28)= 11.791619352031
	sigma_const(13,29)= 13.580888692064
	sigma_const(13,30)= 15.089628980805
	sigma_const(14,16)=  6.450893362682
	sigma_const(14,17)=  8.937436405845
	sigma_const(14,18)= 10.854943255376
	sigma_const(14,19)= 12.052581379413
	sigma_const(14,20)= 12.531951665728
	sigma_const(14,21)= 12.079555422719
	sigma_const(14,22)= 10.434228109502
	sigma_const(14,23)=  8.128995436233
	sigma_const(14,24)=  7.349979562758
	sigma_const(14,25)=  6.381477540041
	sigma_const(14,26)=  6.302120976769
	sigma_const(14,27)=  7.102876211856
	sigma_const(14,28)=  9.390041916469
	sigma_const(14,29)= 11.628597952252
	sigma_const(14,30)= 13.740648056745
	sigma_const(15,17)=  6.353859189371
	sigma_const(15,18)=  8.966864764167
	sigma_const(15,19)= 10.977643424635
	sigma_const(15,20)= 12.340493707655
	sigma_const(15,21)= 12.787309828748
	sigma_const(15,22)= 11.975661489188
	sigma_const(15,23)= 10.402747317438
	sigma_const(15,24)=  9.755843061391
	sigma_const(15,25)=  8.408535806601
	sigma_const(15,26)=  7.054152973254
	sigma_const(15,27)=  6.155958487896
	sigma_const(15,28)=  7.655143121469
	sigma_const(15,29)=  9.837153676797
	sigma_const(15,30)= 12.295932383248
	sigma_const(16,18)=  6.455060083515
	sigma_const(16,19)=  9.100958274589
	sigma_const(16,20)= 11.237797470964
	sigma_const(16,21)= 12.538970140123
	sigma_const(16,22)= 12.620474508707
	sigma_const(16,23)= 11.884013927159
	sigma_const(16,24)= 11.480832843553
	sigma_const(16,25)= 10.220261730983
	sigma_const(16,26)=  8.375731865223
	sigma_const(16,27)=  6.243728825605
	sigma_const(16,28)=  6.120619354446
	sigma_const(16,29)=  7.585913335817
	sigma_const(16,30)= 10.053570657551
	sigma_const(17,19)=  6.416234080267
	sigma_const(17,20)=  9.142999719424
	sigma_const(17,21)= 11.196049306464
	sigma_const(17,22)= 12.139307178662
	sigma_const(17,23)= 12.316574761522
	sigma_const(17,24)= 12.450803755918
	sigma_const(17,25)= 11.756260263942
	sigma_const(17,26)= 10.197794174024
	sigma_const(17,27)=  7.944977942782
	sigma_const(17,28)=  6.432725285681
	sigma_const(17,29)=  6.238905738329
	sigma_const(17,30)=  7.888018780737
	sigma_const(18,20)=  6.459665480289
	sigma_const(18,21)=  9.127011245905
	sigma_const(18,22)= 10.868819875995
	sigma_const(18,23)= 11.948109464512
	sigma_const(18,24)= 12.639189921564
	sigma_const(18,25)= 12.644888039052
	sigma_const(18,26)= 11.667168204526
	sigma_const(18,27)=  9.836439947097
	sigma_const(18,28)=  7.872688075011
	sigma_const(18,29)=  6.303639169085
	sigma_const(18,30)=  6.367968477755
	sigma_const(19,21)=  6.435108776929
	sigma_const(19,22)=  8.880534077633
	sigma_const(19,23)= 10.776194758446
	sigma_const(19,24)= 11.893926014774
	sigma_const(19,25)= 12.599988443633
	sigma_const(19,26)= 12.330401004255
	sigma_const(19,27)= 11.184343587063
	sigma_const(19,28)=  9.363427255114
	sigma_const(19,29)=  7.333791061873
	sigma_const(19,30)=  5.989835075559
	sigma_const(20,22)=  6.348408920539
	sigma_const(20,23)=  8.961451474575
	sigma_const(20,24)= 10.570701901795
	sigma_const(20,25)= 12.065038895019
	sigma_const(20,26)= 12.652689973638
	sigma_const(20,27)= 12.368080805140
	sigma_const(20,28)= 11.131983849897
	sigma_const(20,29)=  9.412403534699
	sigma_const(20,30)=  7.731104416102
	sigma_const(21,23)=  6.473904676176
	sigma_const(21,24)=  8.482780519962
	sigma_const(21,25)= 10.711989848583
	sigma_const(21,26)= 12.160868101194
	sigma_const(21,27)= 12.778702867800
	sigma_const(21,28)= 12.337657191624
	sigma_const(21,29)= 11.290076739868
	sigma_const(21,30)=  9.983538381594
	sigma_const(22,24)=  5.961586665601
	sigma_const(22,25)=  8.819909604302
	sigma_const(22,26)= 11.004507337044
	sigma_const(22,27)= 12.398555634063
	sigma_const(22,28)= 12.855406683461
	sigma_const(22,29)= 12.674003548005
	sigma_const(22,30)= 12.122090344828
	sigma_const(23,25)=  6.480724435539
	sigma_const(23,26)=  9.149680666035
	sigma_const(23,27)= 11.158900999078
	sigma_const(23,28)= 12.477036864953
	sigma_const(23,29)= 13.179002646469
	sigma_const(23,30)= 13.456128178856
	sigma_const(24,26)=  6.475525686900
	sigma_const(24,27)=  9.108983374752
	sigma_const(24,28)= 10.982468045362
	sigma_const(24,29)= 12.289408107097
	sigma_const(24,30)= 13.100425371318
	sigma_const(25,27)=  6.493204685140
	sigma_const(25,28)=  8.989830602548
	sigma_const(25,29)= 11.035878619721
	sigma_const(25,30)= 12.593682786123
	sigma_const(26,28)=  6.329278011430
	sigma_const(26,29)=  9.014680664363
	sigma_const(26,30)= 11.284939706855
	sigma_const(27,29)=  6.567283159008
	sigma_const(27,30)=  9.460645965051
	sigma_const(28,30)=  6.592429526865
	
	
	eps_const(1,1)=40.0
	eps_const(1,2)=30.0
	eps_const(1,3)=20.0
	eps_const(1,4)=17.0
	eps_const(2,1)=30.0
	eps_const(2,2)=25.0
	eps_const(2,3)=13.0
	eps_const(2,4)=10.0
	eps_const(3,1)=20.0
	eps_const(3,2)=13.0
	eps_const(3,3)=5.0
	eps_const(3,4)=2.0
	eps_const(4,1)=17.0
	eps_const(4,2)=10.0
	eps_const(4,3)=2.0
	eps_const(4,4)=1.0
	ENDIF
	

RETURN
END SUBROUTINE SECUENCIAS
