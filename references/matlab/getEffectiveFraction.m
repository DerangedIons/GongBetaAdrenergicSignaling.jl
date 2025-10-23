function [fICaLP,fIKsP,fPLBP,fTnIP,fINaP,fINaKP,fRyRP,fIKurP] = getEffectiveFraction(y, c)
%
% adapted from:  https://github.com/JQXGong/Ohara-beta-adrenergic
%
% Quantitative analysis of variability in an integrated model of
% human ventricular electrophysiology and B-adrenergic signaling
% by Jingqi Q.X. Gong, Monica E. Susilo, Anna Sher, Cynthia J. Musante, Eric A. Sobie
% DOI:https://doi.org/10.1016/j.yjmcc.2020.04.009
%

% Calculating effective fraction of phosphorylated substrates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ICaL
ICaLp = y(40);
fp_ICaL = (ICaLp + c(163)) / c(156); % Fraction of phosphorylated ICaL channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_ICaL = ifthenelse((fp_ICaL < 0.0) ,0.0001 , fp_ICaL);
fp_ICaL = ifthenelse((fp_ICaL > 1.0) ,0.9999 , fp_ICaL);

ical_f_hat_val = (fp_ICaL - c(166)) / (0.9273 - c(166));
ical_f_hat_val = ifthenelse((ical_f_hat_val < 0.0) , 0.0 , ical_f_hat_val); % Effective fraction of phosphorylated ICaL channels
ical_f_hat = ifthenelse((ical_f_hat_val > 1.0) , 1.0 , ical_f_hat_val); % Effective fraction of phosphorylated ICaL channels
%% IKs
IKsp = y(41);
fp_iks = (IKsp + c(145)) / c(144);% Fraction of phosphorylated IKs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_iks = ifthenelse((fp_iks < 0.0) ,0.0001 , fp_iks);
fp_iks = ifthenelse((fp_iks > 1.0) ,0.9999 , fp_iks);

iks_f_hat_val = (fp_iks - c(167)) / (0.785 - c(167));
iks_f_hat_val = ifthenelse((iks_f_hat_val < 0.0) , 0.0 , iks_f_hat_val);% Effective fraction of phosphorylated IKs channels
iks_f_hat = ifthenelse((iks_f_hat_val > 1.0) , 1.0 , iks_f_hat_val);% Effective fraction of phosphorylated IKs channels
%% Iup (PLB)
iup_f_plb = y(42);
% iup_f_pka_val = (iup_f_plb - 0.6591) / (0.9945 - 0.6591);
iup_f_pka_val = (iup_f_plb - 0.6662) / (0.9945 - 0.6662);
iup_f_pka_val = ifthenelse((iup_f_pka_val < 0.0) ,0.0 , iup_f_pka_val);
iup_f_pka = ifthenelse((iup_f_pka_val > 1.0) ,1.0 , iup_f_pka_val);
%% Tni
f_tni = y(43);
calcium_fhat_val = (f_tni - 0.67352) / (0.99918 - 0.67352);
calcium_fhat_val = ifthenelse((calcium_fhat_val < 0.0) , 0.0 , calcium_fhat_val);% Effective fraction of phosphorylated Troponin
calcium_fhat = ifthenelse((calcium_fhat_val > 1.0) , 1.0 , calcium_fhat_val);% Effective fraction of phosphorylated Troponin
%% INa
ina_f_ina = y(44);
ina_f_pka_val = (ina_f_ina - 0.23948) / (0.95014 - 0.23948);
% % ina_f_pka = ifthenelse((ina_f_pka_val < 0.0) , 0.0 , ina_f_pka_val);% Effective fraction of phosphorylated INa channels
ina_f_pka_val = ifthenelse((ina_f_pka_val < 0.0) , 0.0 , ina_f_pka_val);% Effective fraction of phosphorylated INa channels
ina_f_pka = ifthenelse((ina_f_pka_val > 1.0) , 1.0 , ina_f_pka_val);% Effective fraction of phosphorylated INa channels
%% INaK
f_inak = y(45);
inak_fhat_val = (f_inak - 0.12635) / (0.99801 - 0.12635);
inak_fhat_val = ifthenelse((inak_fhat_val < 0.0) , 0.0 , inak_fhat_val);% Effective fraction of phosphorylated INaK pumps
inak_fhat = ifthenelse((inak_fhat_val > 1.0) , 1.0 , inak_fhat_val);% Effective fraction of phosphorylated INaK pumps
%% RyR
RyRp = y(46);
fp_RyR = (RyRp + c(161)) / c(151);% Fraction of phosphorylated RyR channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_RyR = ifthenelse((fp_RyR < 0.0) ,0.0001 , fp_RyR);
fp_RyR = ifthenelse((fp_RyR > 1.0) ,0.9999 , fp_RyR);

irel_fhat_val = (fp_RyR - c(165)) / (0.9586 - c(165));
irel_fhat_val = ifthenelse((irel_fhat_val < 0.0) , 0.0 , irel_fhat_val);% Effective fraction of phosphorylated ryr channels
irel_fhat = ifthenelse((irel_fhat_val > 1.0) , 1.0 , irel_fhat_val);% Effective fraction of phosphorylated ryr channels
%% IKur
f_ikur = y(47);
ikur_fhat_val = (f_ikur -  5.89380e-02) / (0.39375 -  5.89380e-02);
ikur_fhat_val = ifthenelse((ikur_fhat_val < 0.0) , 0.0 , ikur_fhat_val);% Effective fraction of phosphorylated IKur channels
ikur_fhat = ifthenelse((ikur_fhat_val > 1.0) , 1.0 , ikur_fhat_val);% Effective fraction of phosphorylated IKur channels

% Applying the effective fraction of phosphorylated substrates into
% settings struct...  Note: the commented equations below are from the
% code downloaded from Rudy website... The ones above are from Myokit

fICaLP = ical_f_hat;%min(max(((y(98) + sigdata.ICaL_AKAP_PKA) / sigdata.ICaL_tot - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)) / (0.9273 - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)), 0), 1);
fIKsP = iks_f_hat;%min(max(((y(99) + sigdata.IKs_AKAP_PKA) / sigdata.IKs_tot - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)) / (0.7850 - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)), 0), 1);
fPLBP = iup_f_pka;%min(max((y(100) - 6.59100e-001) / (9.94500e-001 - 6.59100e-001), 0), 1);
fTnIP = calcium_fhat;%min(max((y(101) - 6.73519e-001) / (9.99180e-001 - 6.73519e-001), 0), 1);
fINaP = ina_f_pka;%min(max((y(102) - 2.39480e-001) / (9.50143e-001 - 2.39480e-001), 0), 1);
fINaKP = inak_fhat;%min(max((y(103) - 1.26345e-001) / (9.98014e-001 - 1.26345e-001), 0), 1);
fRyRP = irel_fhat;%min(max(((y(104) + sigdata.RyR_AKAP_PKA) / sigdata.RyR_tot - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)) / (0.9586 - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)), 0), 1);
fIKurP = ikur_fhat;%min(max((y(105) - 5.89380e-002) / (3.93747e-001 - 5.89380e-002), 0), 1);

end