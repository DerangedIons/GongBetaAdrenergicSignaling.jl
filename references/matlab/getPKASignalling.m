function ydot = getPKASignalling(y, c)  % Remove t and pace because they are not used....

% adapted from:  https://github.com/JQXGong/Ohara-beta-adrenergic
%
% Quantitative analysis of variability in an integrated model of
% human ventricular electrophysiology and B-adrenergic signaling
% by Jingqi Q.X. Gong, Monica E. Susilo, Anna Sher, Cynthia J. Musante, Eric A. Sobie
% DOI:https://doi.org/10.1016/j.yjmcc.2020.04.009
%

% Create derivatives vector
ydot = zeros(size(y,1), size(y,2));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of the Signaling state variables
beta_cav_Gs_aGTP = y(1);   % 1: Gs_aGTP_CAV
beta_eca_Gs_aGTP = y(2);   % 2: Gs_aGTP_ECAV
beta_cyt_Gs_aGTP = y(3);   % 3: Gs_a_GTP_CYT
beta_cav_Gs_bg = y(4);     % 4: Gs_bg_CAV
beta_eca_Gs_bg = y(5);     % 5: Gs_bg_ECAV
beta_cyt_Gs_bg = y(6);     % 6: Gs_bg_CYT
beta_cav_Gs_aGDP = y(7);   % 7: Gs_aGDP_CAV
beta_eca_Gs_aGDP = y(8);   % 8: Gs_aGDP_ECAV
beta_cyt_Gs_aGDP = y(9);   % 9: Gs_aGDP_CYT

cAMP_cav = y(10);   % 10: cAMP_CAVVV
cAMP_eca = y(11);   % 11: cAMP_ECAV
cAMP_cyt = y(12);  % 12: cAMP_CYT

beta_cav_Rb1_pka_tot = y(13);  % 13: R_pkap_tot_CAV
beta_eca_Rb1_pka_tot = y(14);  % 14: R_pkap_tot_ECAV
beta_cyt_Rb1_pka_tot = y(15);  % 15: R_pkap_tot_CYT
beta_cav_Rb1_grk_tot = y(16);  % 16: R_grkp_tot_CAV
beta_eca_Rb1_grk_tot = y(17);  % 17: R_grkp_tot_ECAV
beta_cyt_Rb1_grk_tot = y(18);  % 18: R_grkp_tot_CYT

pka_cav_ARC = y(19);   % 19: RLC_CAV
pka_cav_A2RC = y(20);  % 20: L2RC_CAV
pka_cav_A2R = y(21);   % 21: L2R_CAV
pka_cav_C = y(22);     % 22: C_CAV
pka_cav_PKIC = y(23);  % 23: PKI_CAV
pka_eca_ARC = y(24);   % 24: RLC_ECAV
pka_eca_A2RC = y(25);  % 25: L2RC_ECAV
pka_eca_A2R = y(26);   % 26: L2R_ECAV
pka_eca_C = y(27);     % 27: C_ECAV
pka_eca_PKIC = y(28);  % 28: PKI_ECAV
pka_cyt_ARC = y(29);   % 29: RLC_CYT
pka_cyt_A2RC = y(30);  % 30: L2RC_CYT
pka_cyt_A2R = y(31);   % 31: L2R_CYT
pka_cyt_C = y(32);     % 32: C_CYT
pka_cyt_PKIC = y(33);  % 33: PKI_CYT

PDE3_P_cav = y(34);    %34   34: PDE3_P_CAV
PDE3_P_cyt = y(35);    %35   35: PDE3_P_CYT
PDE4_P_cav = y(36);    %36   36: PDE4_P_CAV
PDE4_P_eca = y(37);    %37   37: PDE4_P_ECAV
PDE4_P_cyt = y(38);    %38   38: PDE4_P_CYT

inhib1_p = y(39);  %39       39: Inhib1_P_CYT

ICaLp = y(40);     % 40: fLCC_P
IKsp = y(41);      % 41: fIKS_P
iup_f_plb = y(42); % 42: fPLB_P
f_tni = y(43);     % 43: fTnI_P
ina_f_ina = y(44); % 44: fINa_P
f_inak = y(45);    % 45: fINaK_P
RyRp = y(46);      % 46: fRyR_P
f_ikur = y(47);    % 47: fIKur_P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% update here to make the original PKA fractions bound between [0,1]
%%%% use 0.0001 and 0.9999
%%%% for AKAP related 3 channels implement the constraint in EffectiveFraction function
%%% ICaLp IKsp RyRp

%%% iup_f_plb
if iup_f_plb < 0.0
    iup_f_plb=0.0001;
elseif iup_f_plb > 1.0
    iup_f_plb= 0.9999;
end

%%% f_tni
if f_tni  < 0.0
    f_tni =0.0001;
elseif f_tni  > 1.0
    f_tni = 0.9999;
end

%%% ina_f_ina
if ina_f_ina  < 0.0
    ina_f_ina =0.0001;
elseif ina_f_ina > 1.0
    ina_f_ina = 0.9999;
end

%%% f_inak
if f_inak  < 0.0
    f_inak =0.0001;
elseif f_inak > 1.0
    f_inak = 0.9999;
end

%%% f_ikur
if f_ikur  < 0.0
    f_ikur =0.0001;
elseif f_ikur > 1.0
    f_ikur= 0.9999;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta_cav_Rb2_pka_tot = y(48);  % 48: Rb2_pkap_tot_CAV
beta_cav_Rb2_grk_tot = y(49);  % 49: Rb2_grkp_tot_CAV
beta_cav_Gi_aGTP = y(50);      % 50: Gi_aGTP_CAV
beta_cav_Gi_bg = y(51);        % 51: Gi_bg_CAV
beta_cav_Gi_aGDP = y(52);      % 52: Gi_aGDP_CAV
beta_eca_Rb2_pka_tot = y(53);  % 53: Rb2_pkap_tot_ECAV
beta_eca_Rb2_grk_tot = y(54);  % 54: Rb2_grkp_tot_ECAV
beta_eca_Gi_aGTP = y(55);      % 55: Gi_aGTP_ECAV
beta_eca_Gi_bg = y(56);        % 56: Gi_bg_ECAV
beta_eca_Gi_aGDP = y(57);      % 57: Gi_aGDP_ECAV

%%

%% CAVEOLAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total concentration of non-phosphorylated B1AR in the caveolar subspace
beta_cav_Rb1_np_tot = c(87) - beta_cav_Rb1_pka_tot - beta_cav_Rb1_grk_tot;
% Total concentration of non-phosphorylated B2AR in the caveolar subspace
beta_cav_Rb2_np_tot = c(89) - beta_cav_Rb2_pka_tot - beta_cav_Rb2_grk_tot;
% Concentration of Gi holoenzyme in the caveolar subspace
beta_cav_Gi_abg = c(63) * c(59) * c(5) - beta_cav_Gi_aGTP - beta_cav_Gi_aGDP;
% Concentration of Gs holoenzyme in the caveolar subspace
beta_cav_Gs_abg = c(61) * c(58) * c(5) - beta_cav_Gs_aGTP - beta_cav_Gs_aGDP;

beta_cav_Gs_f_d = beta_cav_Gs_abg * c(93) / c(92);
beta_cav_Gs_f_b = (c(94) + c(91)) / c(92) + beta_cav_Rb1_np_tot + beta_cav_Rb2_np_tot - beta_cav_Gs_abg;
beta_cav_Gs_f_c = (c(91) * (beta_cav_Rb1_np_tot - beta_cav_Gs_abg) + c(94) * (beta_cav_Rb2_np_tot - beta_cav_Gs_abg) + c(93)) / c(92);
beta_cav_Gs_f_rr = -beta_cav_Gs_f_d / 27.0 * beta_cav_Gs_f_b ^ 3.0 - beta_cav_Gs_f_b * beta_cav_Gs_f_b * beta_cav_Gs_f_c * beta_cav_Gs_f_c / 108.0 + beta_cav_Gs_f_b * beta_cav_Gs_f_c * beta_cav_Gs_f_d / 6.0 + beta_cav_Gs_f_c ^ 3.0 / 27.0 + beta_cav_Gs_f_d * beta_cav_Gs_f_d / 4.0;
beta_cav_Gs_f_yr = ifthenelse((beta_cav_Gs_f_rr > 0.0) , sqrt(beta_cav_Gs_f_rr) , 0.0) + beta_cav_Gs_f_d / 2.0 + beta_cav_Gs_f_b * beta_cav_Gs_f_c / 6.0 - beta_cav_Gs_f_b ^ 3.0 / 27.0;
beta_cav_Gs_f_yi = ifthenelse((beta_cav_Gs_f_rr < 0.0) , sqrt(-beta_cav_Gs_f_rr) , 0.0);
beta_cav_Gs_f_mag = (beta_cav_Gs_f_yr * beta_cav_Gs_f_yr + beta_cav_Gs_f_yi * beta_cav_Gs_f_yi) ^ (1.0 / 6.0);
beta_cav_Gs_f_arg = atan(beta_cav_Gs_f_yi / beta_cav_Gs_f_yr) / 3.0;
beta_cav_Gs_f_x = (beta_cav_Gs_f_c / 3.0 - beta_cav_Gs_f_b * beta_cav_Gs_f_b / 9.0) / (beta_cav_Gs_f_mag * beta_cav_Gs_f_mag);
beta_cav_Gs_f_r = beta_cav_Gs_f_mag * cos(beta_cav_Gs_f_arg) * (1.0 - beta_cav_Gs_f_x) - beta_cav_Gs_f_b / 3.0;
beta_cav_Gs_f_i = beta_cav_Gs_f_mag * sin(beta_cav_Gs_f_arg) * (1.0 + beta_cav_Gs_f_x);

% Concentration of free Gs in the caveolar subspace
beta_cav_Gs_f = sqrt(beta_cav_Gs_f_r * beta_cav_Gs_f_r + beta_cav_Gs_f_i * beta_cav_Gs_f_i);
% Concentration of free non-phosphorylated beta1AR in caveolar subspace
beta_cav_Rb1_f = beta_cav_Rb1_np_tot / (1.0 + c(1) / c(65) + beta_cav_Gs_f * (c(67) + c(1)) / (c(66) * c(67)));
% Concentration of non-phosphorylated Ligand / Receptor1 complexes in caveolar subspace
beta_cav_LRb1 = c(1) * beta_cav_Rb1_f / c(65);
% Concentration of non-phosphorylated Ligand / Receptor1 / G-protein complexes in caveolar subspace
beta_cav_LRb1Gs = c(1) * beta_cav_Rb1_f * beta_cav_Gs_f / (c(66) * c(67));
% Concentration of free non-phosphorylated beta2AR in caveolar subspace
beta_cav_Rb2_f = beta_cav_Rb2_np_tot / (1.0 + c(1) / c(72) + beta_cav_Gs_f * (c(69) + c(1)) / (c(71) * c(69)));
% Concentration of non-phosphorylated Ligand / Receptor2 complexes in caveolar subspace
beta_cav_LRb2 = c(1) * beta_cav_Rb2_f / c(72);

% Concentration of non-phosphorylated Receptor2 / G-protein complexes in caveolar subspace
beta_cav_Rb2Gs = beta_cav_Rb2_f * beta_cav_Gs_f / c(71);
% Concentration of non-phosphorylated Receptor1 / G-protein complexes in caveolar subspace
beta_cav_Rb1Gs = beta_cav_Rb1_f * beta_cav_Gs_f / c(66);
% Concentration of non-phosphorylated Ligand / Receptor2 / G-protein complexes in caveolar subspace
beta_cav_LRb2Gs = c(1) * beta_cav_Rb2_f * beta_cav_Gs_f / (c(71) * c(69));

% Concentration of total PKA-phosphorylated beta1 receptors
ydot(13) = 0.001 * (c(84) * pka_cav_C * beta_cav_Rb1_np_tot - c(85) * beta_cav_Rb1_pka_tot);
% Concentration of total GRK-phosphorylated beta1 receptors
ydot(16) = 0.001 * (c(83) * c(86) * (beta_cav_LRb1 + beta_cav_LRb1Gs) - c(82) * beta_cav_Rb1_grk_tot);
% Concentration of total PKA-phosphorylated beta2 receptors
ydot(48) = 0.001 * (c(84) * pka_cav_C * beta_cav_Rb2_np_tot - c(85) * beta_cav_Rb2_pka_tot);
% Concentration of total GRK-phosphorylated beta2 receptors
ydot(49) = 0.001 * (c(83) * c(86) * (beta_cav_LRb2 + beta_cav_LRb2Gs) - c(82) * beta_cav_Rb2_grk_tot);

%% EXTRACAVEOLAR %%%%%%%%%%%%%%%%%%
%
% beta_eca
%
% Concentration of Gs holoenzyme in the extracaveolar space
beta_eca_Gs_abg = c(60) * c(58) * c(6) - beta_eca_Gs_aGTP - beta_eca_Gs_aGDP;
% Total concentration of non-phosphorylated B2AR in the extracaveolar space
beta_eca_Rb2_np_tot = c(97) - beta_eca_Rb2_pka_tot - beta_eca_Rb2_grk_tot;
% Total concentration of non-phosphorylated B1AR in the extracaveolar space
beta_eca_Rb1_np_tot = c(98) - beta_eca_Rb1_pka_tot - beta_eca_Rb1_grk_tot;
beta_eca_Gs_f_d = beta_eca_Gs_abg * c(101) / c(102);
beta_eca_Gs_f_c = (c(103) * (beta_eca_Rb1_np_tot - beta_eca_Gs_abg) + c(100) * (beta_eca_Rb2_np_tot - beta_eca_Gs_abg) + c(101)) / c(102);
beta_eca_Gs_f_b = (c(100) + c(103)) / c(102) + beta_eca_Rb1_np_tot + beta_eca_Rb2_np_tot - beta_eca_Gs_abg;
beta_eca_Gs_f_rr = -beta_eca_Gs_f_d / 27.0 * beta_eca_Gs_f_b ^ 3.0 - beta_eca_Gs_f_b * beta_eca_Gs_f_b * beta_eca_Gs_f_c * beta_eca_Gs_f_c / 108.0 + beta_eca_Gs_f_b * beta_eca_Gs_f_c * beta_eca_Gs_f_d / 6.0 + beta_eca_Gs_f_c ^ 3.0 / 27.0 + beta_eca_Gs_f_d * beta_eca_Gs_f_d / 4.0;
beta_eca_Gs_f_yi = ifthenelse((beta_eca_Gs_f_rr < 0.0) , sqrt(-beta_eca_Gs_f_rr) , 0.0);
beta_eca_Gs_f_yr = ifthenelse((beta_eca_Gs_f_rr > 0.0) , sqrt(beta_eca_Gs_f_rr) , 0.0) + beta_eca_Gs_f_d / 2.0 + beta_eca_Gs_f_b * beta_eca_Gs_f_c / 6.0 - beta_eca_Gs_f_b ^ 3.0 / 27.0;
beta_eca_Gs_f_mag = (beta_eca_Gs_f_yr * beta_eca_Gs_f_yr + beta_eca_Gs_f_yi * beta_eca_Gs_f_yi) ^ (1.0 / 6.0);
beta_eca_Gs_f_arg = atan(beta_eca_Gs_f_yi / beta_eca_Gs_f_yr) / 3.0;
beta_eca_Gs_f_x = (beta_eca_Gs_f_c / 3.0 - beta_eca_Gs_f_b * beta_eca_Gs_f_b / 9.0) / (beta_eca_Gs_f_mag * beta_eca_Gs_f_mag);
beta_eca_Gs_f_i = beta_eca_Gs_f_mag * sin(beta_eca_Gs_f_arg) * (1.0 + beta_eca_Gs_f_x);
beta_eca_Gs_f_r = beta_eca_Gs_f_mag * cos(beta_eca_Gs_f_arg) * (1.0 - beta_eca_Gs_f_x) - beta_eca_Gs_f_b / 3.0;
% Concentration of free Gs in the caveolar subspace
beta_eca_Gs_f = sqrt(beta_eca_Gs_f_r * beta_eca_Gs_f_r + beta_eca_Gs_f_i * beta_eca_Gs_f_i);
% Concentration of free non-phosphorylated beta1AR in the extracaveolar space
beta_eca_Rb1_f = beta_eca_Rb1_np_tot / (1.0 + c(1) / c(65) + beta_eca_Gs_f * (c(67) + c(1)) / (c(66) * c(67)));
% Concentration of free non-phosphorylated beta2AR in the extracaveolar space
beta_eca_Rb2_f = beta_eca_Rb2_np_tot / (1.0 + c(1) / c(72) + beta_eca_Gs_f * (c(69) + c(1)) / (c(71) * c(69)));
% Concentration of non-phosphorylated Ligand / Receptor1 complexes in the extracaveolar space
beta_eca_LRb1 = c(1) * beta_eca_Rb1_f / c(65);
% Concentration of non-phosphorylated Ligand / Receptor2 complexes in the extracaveolar space
beta_eca_LRb2 = c(1) * beta_eca_Rb2_f / c(72);
% Concentration of non-phosphorylated Ligand / Receptor2 / G-protein complexes in the extracaveolar space
beta_eca_LRb2Gs = c(1) * beta_eca_Rb2_f * beta_eca_Gs_f / (c(71) * c(69));
% Concentration of non-phosphorylated Ligand / Receptor1 / G-protein complexes in the extracaveolar space
beta_eca_LRb1Gs = c(1) * beta_eca_Rb1_f * beta_eca_Gs_f / (c(66) * c(67));
% Concentration of non-phosphorylated Receptor2 / G-protein complexes in the extracaveolar space
beta_eca_Rb2Gs = beta_eca_Rb2_f * beta_eca_Gs_f / c(71);
% Concentration of non-phosphorylated Receptor1 / G-protein complexes in the extracaveolar space
beta_eca_Rb1Gs = beta_eca_Rb1_f * beta_eca_Gs_f / c(66);
beta_eca_RGs_tot = beta_eca_Rb1Gs + c(96) * beta_eca_Rb2Gs;
beta_eca_LRGs_tot = beta_eca_LRb1Gs + c(96) * beta_eca_LRb2Gs;
% Concentration of Gi holoenzyme in the extracaveolar space
beta_eca_Gi_abg = c(64) * c(59) * c(6) - beta_eca_Gi_aGTP - beta_eca_Gi_aGDP;


beta_eca_Rb2_pka_f_c = -beta_eca_Rb2_pka_tot * c(70) * c(73);
beta_eca_Rb2_pka_f_b = beta_eca_Gi_abg * (c(1) + c(73)) - beta_eca_Rb2_pka_tot * (c(73) + c(1)) + c(70) * c(73) * (1.0 + c(1) / c(68));
% Concentration of free PKA-phosphorylated beta2AR in the extracaveolar space
beta_eca_Rb2_pka_f = (-beta_eca_Rb2_pka_f_b + sqrt(beta_eca_Rb2_pka_f_b * beta_eca_Rb2_pka_f_b - 4.0 * c(99) * beta_eca_Rb2_pka_f_c)) / (2.0 * c(99));

% Concentration of total PKA-phosphorylated beta1 receptors in the extracaveolar space
ydot(14) = 0.001 * (c(84) * pka_eca_C * beta_eca_Rb1_np_tot - c(85) * beta_eca_Rb1_pka_tot);
% Concentration of total GRK-phosphorylated beta1 receptors in the extracaveolar space
ydot(17) = 0.001 * (c(83) * c(95) * (beta_eca_LRb1 + beta_eca_LRb1Gs) - c(82) * beta_eca_Rb1_grk_tot);
% Concentration of total PKA-phosphorylated beta2 receptors in the extracaveolar space
ydot(53) = 0.001 * (c(84) * pka_eca_C * beta_eca_Rb2_np_tot - c(85) * beta_eca_Rb2_pka_tot);
% Concentration of total GRK-phosphorylated beta2 receptors in the extracaveolar space
ydot(54) = 0.001 * (c(83) * c(95) * (beta_eca_LRb2 + beta_eca_LRb2Gs) - c(82) * beta_eca_Rb2_grk_tot);


%% CYTOPLASM %%%%%%%%%%%
% Concentration of Gs holoenzyme in the cytoplasm
beta_cyt_Gs_abg = c(62) * c(58) * c(7) - beta_cyt_Gs_aGTP - beta_cyt_Gs_aGDP;
% Total concentration of non-phosphorylated beta-1 AR in the cytoplasm
beta_cyt_Rb1_np_tot = c(104) - beta_cyt_Rb1_pka_tot - beta_cyt_Rb1_grk_tot;

beta_cyt_Rb1_np_f_b = beta_cyt_Gs_abg * (c(67) + c(1)) - beta_cyt_Rb1_np_tot * (c(67) + c(1)) + c(66) * c(67) * (1.0 + c(1) / c(65));
beta_cyt_Rb1_np_f_c = -beta_cyt_Rb1_np_tot * c(67) * c(66);
% Concentration of free non-phosphorylated beta-1 AR in the cytoplasm
Rb1_np_f = (-beta_cyt_Rb1_np_f_b + sqrt(beta_cyt_Rb1_np_f_b * beta_cyt_Rb1_np_f_b - 4.0 * c(106) * beta_cyt_Rb1_np_f_c)) / (2.0 * c(106));
% Concentration of free (non-complexed) Gi in the cytoplasm
beta_cyt_Gs_f = beta_cyt_Gs_abg / (1.0 + Rb1_np_f / c(66) * (1.0 + c(1) / c(67)));
% Concentration of non-phosphorylated ligand / beta-1 AR complexes in the cytoplasm
LRb1_np = c(1) * Rb1_np_f / c(65);
% Concentration of non-phosphorylated ligand / beta-1 AR / Gs complexes in the cytoplasm
LRb1Gs_np = c(1) * Rb1_np_f * beta_cyt_Gs_f / (c(66) * c(67));
% Concentration of non-phosphorylated beta-1 AR / Gs complexes in the cytoplasm
Rb1Gs_np = beta_cyt_Gs_f * Rb1_np_f / c(66);

% Concentration of total PKA-phosphorylated receptors in the cytoplasm
ydot(15) = 0.001 * (c(84) * pka_cyt_C * beta_cyt_Rb1_np_tot - c(85) * beta_cyt_Rb1_pka_tot);
% Concentration of total GRK-phosphorylated receptors in the cytoplasm
ydot(18) = 0.001 * (c(83) * c(105) * (LRb1_np + LRb1Gs_np) - c(82) * beta_cyt_Rb1_grk_tot);


%% function Mod_GprotAct() in the other version %%%

beta_cav_RGs_tot = beta_cav_Rb1Gs + c(88) * beta_cav_Rb2Gs;
beta_cav_LRGs_tot = beta_cav_LRb1Gs + c(88) * beta_cav_LRb2Gs;
beta_cav_Rb2_pka_f_c = -beta_cav_Rb2_pka_tot * c(70) * c(73);
beta_cav_Rb2_pka_f_b = beta_cav_Gi_abg * (c(1) + c(73)) - beta_cav_Rb2_pka_tot * (c(73) + c(1)) + c(70) * c(73) * (1.0 + c(1) / c(68));
% Concentration of free PKA-phosphorylated beta2AR in the caveolar subspace
beta_cav_Rb2_pka_f = (-beta_cav_Rb2_pka_f_b + sqrt(beta_cav_Rb2_pka_f_b * beta_cav_Rb2_pka_f_b - 4.0 * c(90) * beta_cav_Rb2_pka_f_c)) / (2.0 * c(90));
% Concentration of free (non-complexed) Gi in the caveolar subspace
beta_cav_Gi_f = beta_cav_Gi_abg / (1.0 + beta_cav_Rb2_pka_f / c(70) * (1.0 + c(1) / c(73)));
% Concentration of PKA phosphorylated b2AR/Gi complexes in caveolar subspace
beta_cav_Rb2Gi = beta_cav_Rb2_pka_f * beta_cav_Gi_f / c(70);
% Concentration of PKA phosphorylated Ligand/b2AR/Gi complexes in caveolar subspace
beta_cav_LRb2Gi = beta_cav_Rb2Gi * c(1) / c(73);
% Concentration of free (non-complexed) Gi in the extracaveolar space
beta_eca_Gi_f = beta_eca_Gi_abg / (1.0 + beta_eca_Rb2_pka_f / c(70) * (1.0 + c(1) / c(73)));
% Concentration of PKA phosphorylated b2AR/Gi complexes in extracaveolar space
beta_eca_Rb2Gi = beta_eca_Rb2_pka_f * beta_eca_Gi_f / c(70);
% Concentration of PKA phosphorylated Ligand/b2AR/Gi complexes in extracaveolar space
beta_eca_LRb2Gi = c(1) / c(73) * beta_eca_Rb2Gi;

% Concentration of active Gs alpha subunit in caveolar subspace (tri-phosphate)
ydot(1) = 0.001 * (c(75) * beta_cav_RGs_tot + c(74) * beta_cav_LRGs_tot - c(78) * beta_cav_Gs_aGTP);
% Concentration of active Gi alpha subunit in caveolar subspace
ydot(50) = 0.001 * (c(77) * beta_cav_Rb2Gi + c(76) * beta_cav_LRb2Gi - c(79) * beta_cav_Gi_aGTP);
% Concentration of active Gs alpha subunit in the extracaveolar space (tri-phosphate)
ydot(2) = 0.001 * (c(75) * beta_eca_RGs_tot + c(74) * beta_eca_LRGs_tot - c(78) * beta_eca_Gs_aGTP);
% Concentration of active Gi alpha subunit in the extracaveolar space
ydot(55) = 0.001 * (c(77) * beta_eca_Rb2Gi + c(76) * beta_eca_LRb2Gi - c(79) * beta_eca_Gi_aGTP);
% Concentration of active Gs alpha subunit in cytoplasm
ydot(3) = 0.001 * (c(75) * Rb1Gs_np + c(74) * LRb1Gs_np - c(78) * beta_cyt_Gs_aGTP);

% Concentration of active Gs beta-gamma subunit in caveolar subspace
ydot(4) = 0.001 * (c(75) * beta_cav_RGs_tot + c(74) * beta_cav_LRGs_tot - c(80) * beta_cav_Gs_bg * beta_cav_Gs_aGDP);
% Concentration of active Gi beta-gamma subunit in caveolar subspace
ydot(51) = 0.001 * (c(77) * beta_cav_Rb2Gi + c(76) * beta_cav_LRb2Gi - c(81) * beta_cav_Gi_bg * beta_cav_Gi_aGDP);
% Concentration of active Gs beta-gamma subunit in the extracaveolar space
ydot(5) = 0.001 * (c(75) * beta_eca_RGs_tot + c(74) * beta_eca_LRGs_tot - c(80) * beta_eca_Gs_bg * beta_eca_Gs_aGDP);
% Concentration of active Gi beta-gamma subunit in the extracaveolar space
ydot(56) = 0.001 * (c(77) * beta_eca_Rb2Gi + c(76) * beta_eca_LRb2Gi - c(81) * beta_eca_Gi_bg * beta_eca_Gi_aGDP);
% Concentration of active Gs beta-gamma subunit in cytoplasm
ydot(6) = 0.001 * (c(75) * Rb1Gs_np + c(74) * LRb1Gs_np - c(80) * beta_cyt_Gs_bg * beta_cyt_Gs_aGDP);

% Concentration of inactive Gs alpha subunit in caveolar subspace (di-phosphate)
ydot(7) = 0.001 * (c(78) * beta_cav_Gs_aGTP - c(80) * beta_cav_Gs_bg * beta_cav_Gs_aGDP);
% Concentration of inactive Gi alpha subunit in caveolar subspace
ydot(52) = 0.001 * (c(79) * beta_cav_Gi_aGTP - c(81) * beta_cav_Gi_bg * beta_cav_Gi_aGDP);
% Concentration of inactive Gs alpha subunit in the extracaveolar space (di-phosphate)
ydot(8) = 0.001 * (c(78) * beta_eca_Gs_aGTP - c(80) * beta_eca_Gs_bg * beta_eca_Gs_aGDP);
% Concentration of inactive Gi alpha subunit in the extracaveolar space
ydot(57) = 0.001 * (c(79) * beta_eca_Gi_aGTP - c(81) * beta_eca_Gi_bg * beta_eca_Gi_aGDP);
% Concentration of inactive Gs alpha subunit in cytoplasm
ydot(9) = 0.001 * (c(78) * beta_cyt_Gs_aGTP - c(80) * beta_cyt_Gs_bg * beta_cyt_Gs_aGDP);

%% function Mod_AC()

% only calculating constants, but does not compute derivative of state
% variables

%% function Mod_PKA() in the other version %%%%%%%%%%
% Concentration of free PKA RC subunits in the caveolar compartment
pka_cav_RCf = c(12) - pka_cav_ARC - pka_cav_A2RC - pka_cav_A2R;

% Caveolar concentration of PKA RC dimer with 1 cAMP molecule bound
ydot(19) = 0.001 * (c(19) * pka_cav_RCf * cAMP_cav - c(22) * pka_cav_ARC - c(20) * pka_cav_ARC * cAMP_cav + c(23) * pka_cav_A2RC);
% Caveolar concentration of PKA RC dimer with 2 cAMP molecules bound
ydot(20) = 0.001 * (c(20) * pka_cav_ARC * cAMP_cav - (c(23) + c(21)) * pka_cav_A2RC + c(24) * pka_cav_A2R * pka_cav_C);
% Caveolar concentration of PKA R subunit with 2 cAMP molecules bound
ydot(21) = 0.001 * (c(21) * pka_cav_A2RC - c(24) * pka_cav_A2R * pka_cav_C);
% Caveolar concentration of free PKA catalytic subunit
ydot(22) = 0.001 * (c(21) * pka_cav_A2RC - c(24) * pka_cav_A2R * pka_cav_C + c(18) * pka_cav_PKIC - c(17) * (c(14) - pka_cav_PKIC) * pka_cav_C);
% Caveolar concentration of free PKI inactivated PKA C subunit
ydot(23) = 0.001 * (c(17) * (c(14) - pka_cav_PKIC) * pka_cav_C - c(18) * pka_cav_PKIC);

% Concentration of free PKA RC subunits in the Extracaveolar compartment
pka_eca_RCf = c(11) - pka_eca_ARC - pka_eca_A2RC - pka_eca_A2R;
% Extracaveolar rate of change in free cAMP through binding by PKA
pka_eca_dcAMP = -c(19) * pka_eca_RCf * cAMP_eca + c(25) * pka_eca_ARC - c(20) * pka_eca_ARC * cAMP_eca + c(26) * pka_eca_A2RC;

% Extracaveolar concentration of PKA RC dimer with 1 cAMP molecule bound
ydot(24) = 0.001 * (c(19) * pka_eca_RCf * cAMP_eca - c(25) * pka_eca_ARC - c(20) * pka_eca_ARC * cAMP_eca + c(26) * pka_eca_A2RC);
% Extracaveolar concentration of PKA RC dimer with 2 cAMP molecules bound
ydot(25) = 0.001 * (c(20) * pka_eca_ARC * cAMP_eca - (c(26) + c(21)) * pka_eca_A2RC + c(27) * pka_eca_A2R * pka_eca_C);
% Extracaveolar concentration of PKA R subunit with 2 cAMP molecules bound
ydot(26) = 0.001 * (c(21) * pka_eca_A2RC - c(27) * pka_eca_A2R * pka_eca_C);
% Extracaveolar concentration of free PKA catalytic subunit
ydot(27) = 0.001 * (c(21) * pka_eca_A2RC - c(27) * pka_eca_A2R * pka_eca_C + c(18) * pka_eca_PKIC - c(17) * (c(15) - pka_eca_PKIC) * pka_eca_C);
% Extracaveolar concentration of free PKI inactivated PKA C subunit
ydot(28) = 0.001 * (c(17) * (c(15) - pka_eca_PKIC) * pka_eca_C - c(18) * pka_eca_PKIC);

% Concentration of free PKA RC subunits in the Cytosolic compartment
pka_cyt_RCf = c(13) - pka_cyt_ARC - pka_cyt_A2RC - pka_cyt_A2R;

% Cytosolic concentration of PKA RC dimer with 1 cAMP molecule bound
ydot(29) = 0.001 * (c(19) * pka_cyt_RCf * cAMP_cyt - c(28) * pka_cyt_ARC - c(20) * pka_cyt_ARC * cAMP_cyt + c(29) * pka_cyt_A2RC);
% Cytosolic concentration of PKA RC dimer with 2 cAMP molecules bound
ydot(30) = 0.001 * (c(20) * pka_cyt_ARC * cAMP_cyt - (c(29) + c(21)) * pka_cyt_A2RC + c(30) * pka_cyt_A2R * pka_cyt_C);
% Cytosolic concentration of PKA R subunit with 2 cAMP molecules bound
ydot(31) = 0.001 * (c(21) * pka_cyt_A2RC - c(30) * pka_cyt_A2R * pka_cyt_C);
% Cytosolic concentration of free PKA catalytic subunit
ydot(32) = 0.001 * (c(21) * pka_cyt_A2RC - c(30) * pka_cyt_A2R * pka_cyt_C + c(18) * pka_cyt_PKIC - c(17) * (c(16) - pka_cyt_PKIC) * pka_cyt_C);
% Cytosolic concentration of free PKI inactivated PKA C subunit
ydot(33) = 0.001 * (c(17) * (c(16) - pka_cyt_PKIC) * pka_cyt_C - c(18) * pka_cyt_PKIC);

%% function Mod_cAMP() in the other version


% Caveolar rate of change in free cAMP through binding by PKA
pka_cav_dcAMP = -c(19) * pka_cav_RCf * cAMP_cav + c(22) * pka_cav_ARC - c(20) * pka_cav_ARC * cAMP_cav + c(23) * pka_cav_A2RC;

% Cytosolic rate of change in free cAMP through binding by PKA
pka_cyt_dcAMP = -c(19) * pka_cyt_RCf * cAMP_cyt + c(28) * pka_cyt_ARC - c(20) * pka_cyt_ARC * cAMP_cyt + c(29) * pka_cyt_A2RC;

%PDE
% Rate of cAMP degradation by PDE2 in cytosolic subspace
dcAMP_PDE2_cyt = c(52) * c(41) / (1.0 + c(44) / cAMP_cyt);
% Rate of cAMP degradation by PDE2 in extracaveolar subspace
dcAMP_PDE2_eca = c(51) * c(41) / (1.0 + c(44) / cAMP_eca);
% Rate of cAMP degradation by PDE2 in caveolar subspace
dcAMP_PDE2_cav = c(50) * c(41) / (1.0 + c(44) / cAMP_cav);
% Rate of cAMP degradation by PDE3 in caveolar subspace
dcAMP_PDE3_cav = (c(53) + (c(47) - 1.0) * PDE3_P_cav) * c(42) / (1.0 + c(45) / cAMP_cav);
% Rate of cAMP degradation by PDE4 in cytosolic subspace
dcAMP_PDE4_cyt = (c(57) + (c(47) - 1.0) * PDE4_P_cyt) * c(43) / (1.0 + c(46) / cAMP_cyt);
% Rate of cAMP degradation by PDE4 in extracaveolar subspace
dcAMP_PDE4_eca = (c(56) + (c(47) - 1.0) * PDE4_P_eca) * c(43) / (1.0 + c(46) / cAMP_eca);
% Rate of cAMP degradation by PDE4 in caveolar subspace
dcAMP_PDE4_cav = (c(55) + (c(47) - 1.0) * PDE4_P_cav) * c(43) / (1.0 + c(46) / cAMP_cav);
% Rate of cAMP degradation by PDE3 in cytosolic subspace
dcAMP_PDE3_cyt = (c(54) + (c(47) - 1.0) * PDE3_P_cyt) * c(42) / (1.0 + c(45) / cAMP_cyt);

camp_cAMP_cyt_pde = dcAMP_PDE2_cyt + dcAMP_PDE3_cyt + dcAMP_PDE4_cyt;
camp_cAMP_cyt_j1 = c(9) * (cAMP_cav - cAMP_cyt) / c(4);
camp_cAMP_cyt_j2 = c(10) * (cAMP_eca - cAMP_cyt) / c(4);

camp_cAMP_eca_pde = dcAMP_PDE2_eca + dcAMP_PDE4_eca;
camp_cAMP_eca_j2 = c(10) * (cAMP_eca - cAMP_cyt) / c(3);
camp_cAMP_eca_j1 = c(8) * (cAMP_cav - cAMP_eca) / c(3);

camp_cAMP_cav_pde = dcAMP_PDE2_cav + dcAMP_PDE3_cav + dcAMP_PDE4_cav;
camp_cAMP_cav_j2 = c(9) * (cAMP_cav - cAMP_cyt) / c(2);
camp_cAMP_cav_j1 = c(8) * (cAMP_cav - cAMP_eca) / c(2);
% ac
%
ac_kAC47_cyt_gsa = beta_cyt_Gs_aGTP ^ c(107);
kAC47_cyt = c(116) * (c(114) + ac_kAC47_cyt_gsa / (c(110) + ac_kAC47_cyt_gsa));
ac_kAC56_cav_gsa = beta_cav_Gs_aGTP ^ c(108);
gsi = beta_cav_Gs_aGTP ^ c(109);
kAC56_cav = c(117) * (c(115) + ac_kAC56_cav_gsa / (c(111) + ac_kAC56_cav_gsa)) * (1.0 - (1.0 - c(118) * gsi / (c(113) + gsi)) * beta_cav_Gi_bg / (c(112) + beta_cav_Gi_bg));
ac_kAC47_eca_gsa = beta_eca_Gs_aGTP ^ c(107);
kAC47_eca = c(116) * (c(114) + ac_kAC47_eca_gsa / (c(110) + ac_kAC47_eca_gsa));
ac_kAC56_cyt_gsa = beta_cyt_Gs_aGTP ^ c(108);
kAC56_cyt = c(117) * (c(115) + ac_kAC56_cyt_gsa / (c(111) + ac_kAC56_cyt_gsa));

% Rate of cAMP production by AC type 4/7 in cytoplasm
dcAMP_AC47_cyt = kAC47_cyt * c(119) * c(121);
% Rate of cAMP production by AC type 5/6 in cytoplasm
dcAMP_AC56_cyt = kAC56_cyt * c(123) * c(121);
% Rate of cAMP production by AC type 5/6 in caveolar subspace
dcAMP_AC56_cav = kAC56_cav * c(120) * c(121);
% Rate of cAMP production by AC type 4/7 in extracaveolar subspace
dcAMP_AC47_eca = kAC47_eca * c(122) * c(121);


% Caveolar concentration of cAMP
ydot(10) = 0.001 * (pka_cav_dcAMP + dcAMP_AC56_cav - camp_cAMP_cav_pde - camp_cAMP_cav_j1 - camp_cAMP_cav_j2);
% Extracaveolar concentration of cAMP
ydot(11) = 0.001 * (pka_eca_dcAMP + dcAMP_AC47_eca - camp_cAMP_eca_pde + camp_cAMP_eca_j1 - camp_cAMP_eca_j2);
% Cytosolic concentration of cAMP
ydot(12) = 0.001 * (pka_cyt_dcAMP + dcAMP_AC47_cyt + dcAMP_AC56_cyt - camp_cAMP_cyt_pde + camp_cAMP_cyt_j1 + camp_cAMP_cyt_j2);


%% Mod_PDE_Phosphorylation() function

% Concentration of phosphorylated PDE3 in the caveolar subspace
ydot(34) = 0.001 * (c(48) * pka_cav_C * (c(53) - PDE3_P_cav) - c(49) * PDE3_P_cav);
% Concentration of phosphorylated PDE3 in the cytosolic subspace
ydot(35) = 0.001 * (c(48) * pka_cyt_C * (c(54) - PDE3_P_cyt) - c(49) * PDE3_P_cyt);
% Concentration of phosphorylated PDE4 in the caveolar subspace
ydot(36) = 0.001 * (c(48) * pka_cav_C * (c(55) - PDE4_P_cav) - c(49) * PDE4_P_cav);
% Concentration of phosphorylated PDE4 in the extracaveolar subspace
ydot(37) = 0.001 * (c(48) * pka_eca_C * (c(56) - PDE4_P_eca) - c(49) * PDE4_P_eca);
% Concentration of phosphorylated PDE4 in the cytosolic subspace
ydot(38) = 0.001 * (c(48) * pka_cyt_C * (c(57) - PDE4_P_cyt) - c(49) * PDE4_P_cyt);

%% Mod_PP1_Inhibition()

pp1_PP1f_cyt_sum = c(38) - c(37) + inhib1_p;
% Concentration of uninhibited PP1 in the cytosolic compartment
PP1f_cyt = 0.5 * (sqrt(pp1_PP1f_cyt_sum ^ 2.0 + 4.0 * c(38) * c(37)) - pp1_PP1f_cyt_sum);
di = c(39) - inhib1_p;
% Concentration of phosphorylated PP1 inhibitor 1 (cytoplasmic)
ydot(39) = 0.001 * (c(31) * pka_cyt_C * di / (c(33) + di) - c(32) * c(40) * inhib1_p / (c(34) + inhib1_p));

%% Mod_Channel_Phosphorylation()

% Substrates without AKAP
% Fraction of phosphorylated PLB
ydot(42) = 0.001 * (c(132) * pka_cyt_C * (1.0 - iup_f_plb) / (c(134) + 1.0 - iup_f_plb) - c(133) * PP1f_cyt * iup_f_plb / (c(135) + iup_f_plb));
% Fraction of phosphorylated Troponin
ydot(43) = 0.001 * (c(147) * pka_cyt_C * (1.0 - f_tni) / (c(149) + 1.0 - f_tni) - c(148) * c(40) * f_tni / (c(150) + f_tni));
% Fraction of phosphorylated INa channels
ydot(44) = 0.001 * (c(130) * pka_cav_C * (1.0 - ina_f_ina) / (c(128) + 1.0 - ina_f_ina) - c(131) * c(36) * ina_f_ina / (c(129) + ina_f_ina));
% Fraction of phosphorylated INaK
ydot(45) = 0.001 * (c(124) * pka_cav_C * (1.0 - f_inak) / (c(126) + 1.0 - f_inak) - c(125) * c(36) * f_inak / (c(127) + f_inak));
% Fraction of phosphorylated IKur channels
ydot(47) = 0.001 * (c(136) * pka_eca_C * (1.0 - f_ikur) / (c(138) + 1.0 - f_ikur) - c(137) * c(35) * f_ikur / (c(139) + f_ikur));


%Substrates with AKAP
iks_sig_IKsp_dif = c(146) - IKsp;
% Concentration of phosphorylated IKs channels
ydot(41) = 0.001 * (c(140) * pka_eca_C * iks_sig_IKsp_dif / (c(142) + iks_sig_IKsp_dif) - c(141) * c(35) * IKsp / (c(143) + IKsp));
akap_sig_RyRp_dif = c(162) - RyRp;
ydot(46) = 0.001 * (c(152) * pka_cav_C * akap_sig_RyRp_dif / (c(154) + akap_sig_RyRp_dif) - c(153) * c(36) * RyRp / (c(155) + RyRp));
akap_sig_ICaLp_dif = c(164) - ICaLp;
% Concentration of phosphorylated L-type Calcium channels
ydot(40) = 0.001 * (c(157) * pka_cav_C * akap_sig_ICaLp_dif / (c(159) + akap_sig_ICaLp_dif) - c(158) * c(36) * ICaLp / (c(160) + ICaLp));

end