clear all;

nm      = 10^(-9);

lambda0min    = 380*nm;
lambda0max    = 750*nm;
lambda0points = 74;
lambda01      = lambda0min + (lambda0max-lambda0min)*[0:lambda0points]/lambda0points;

%--------------------------------------------------------------TITANIUM/AIR
load('Ti_Air_Data.mat');
load('Air_Ti_Data.mat');
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_Air_PmaxTEabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, Air_Ti_PmaxTEabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('Ti/Air max P^{TE}_{abs}/P^{TE}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_Air_PmaxTMabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, Air_Ti_PmaxTMabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_Air_PmaxTEscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, Air_Ti_PmaxTEscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_Air_PmaxTMscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, Air_Ti_PmaxTMscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_Air_aoboptTEabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_Air_boloptTEabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, Air_Ti_aoboptTEabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, Air_Ti_boloptTEabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('Ti/Air max P^{TE}_{abs}/P^{TE}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_Air_aoboptTMabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_Air_boloptTMabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, Air_Ti_aoboptTMabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, Air_Ti_boloptTMabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_Air_aoboptTEscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_Air_boloptTEscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, Air_Ti_aoboptTEscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, Air_Ti_boloptTEscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_Air_aoboptTMscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_Air_boloptTMscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, Air_Ti_aoboptTMscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, Air_Ti_boloptTMscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);

%-------------------------------------------------------------TITANIUM/ALSB
load('Ti_AlSb_Data.mat');
load('AlSb_Ti_Data.mat');
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_AlSb_PmaxTEabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, AlSb_Ti_PmaxTEabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('Ti/AlSb max P^{TE}_{abs}/P^{TE}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_AlSb_PmaxTMabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, AlSb_Ti_PmaxTMabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_AlSb_PmaxTEscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, AlSb_Ti_PmaxTEscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_AlSb_PmaxTMscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, AlSb_Ti_PmaxTMscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_AlSb_aoboptTEabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_AlSb_boloptTEabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, AlSb_Ti_aoboptTEabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, AlSb_Ti_boloptTEabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('Ti/AlSb max P^{TE}_{abs}/P^{TE}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_AlSb_aoboptTMabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_AlSb_boloptTMabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, AlSb_Ti_aoboptTMabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, AlSb_Ti_boloptTMabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_AlSb_aoboptTEscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_AlSb_boloptTEscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, AlSb_Ti_aoboptTEscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, AlSb_Ti_boloptTEscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_AlSb_aoboptTMscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_AlSb_boloptTMscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, AlSb_Ti_aoboptTMscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, AlSb_Ti_boloptTMscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);

%--------------------------------------------------------------TITANIUM/ASI
load('Ti_aSi_Data.mat');
load('aSi_Ti_Data.mat');
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_aSi_PmaxTEabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, aSi_Ti_PmaxTEabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('Ti/aSi max P^{TE}_{abs}/P^{TE}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_aSi_PmaxTMabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, aSi_Ti_PmaxTMabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_aSi_PmaxTEscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, aSi_Ti_PmaxTEscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_aSi_PmaxTMscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, aSi_Ti_PmaxTMscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_aSi_aoboptTEabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_aSi_boloptTEabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, aSi_Ti_aoboptTEabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, aSi_Ti_boloptTEabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('Ti/aSi max P^{TE}_{abs}/P^{TE}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_aSi_aoboptTMabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_aSi_boloptTMabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, aSi_Ti_aoboptTMabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, aSi_Ti_boloptTMabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_aSi_aoboptTEscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_aSi_boloptTEscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, aSi_Ti_aoboptTEscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, aSi_Ti_boloptTEscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_aSi_aoboptTMscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_aSi_boloptTMscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, aSi_Ti_aoboptTMscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, aSi_Ti_boloptTMscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);

%--------------------------------------------------------------TITANIUM/GAP
load('Ti_GaP_Data.mat');
load('GaP_Ti_Data.mat');
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_GaP_PmaxTEabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, GaP_Ti_PmaxTEabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('Ti/GaP max P^{TE}_{abs}/P^{TE}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_GaP_PmaxTMabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, GaP_Ti_PmaxTMabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_GaP_PmaxTEscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, GaP_Ti_PmaxTEscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_GaP_PmaxTMscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, GaP_Ti_PmaxTMscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_GaP_aoboptTEabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_GaP_boloptTEabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, GaP_Ti_aoboptTEabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, GaP_Ti_boloptTEabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('Ti/GaP max P^{TE}_{abs}/P^{TE}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_GaP_aoboptTMabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_GaP_boloptTMabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, GaP_Ti_aoboptTMabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, GaP_Ti_boloptTMabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_GaP_aoboptTEscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_GaP_boloptTEscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, GaP_Ti_aoboptTEscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, GaP_Ti_boloptTEscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_GaP_aoboptTMscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_GaP_boloptTMscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, GaP_Ti_aoboptTMscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, GaP_Ti_boloptTMscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);

%-------------------------------------------------------------TITANIUM/GASB
load('Ti_GaSb_Data.mat');
load('GaSb_Ti_Data.mat');
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_GaSb_PmaxTEabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, GaSb_Ti_PmaxTEabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('Ti/GaSb max P^{TE}_{abs}/P^{TE}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_GaSb_PmaxTMabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, GaSb_Ti_PmaxTMabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_GaSb_PmaxTEscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, GaSb_Ti_PmaxTEscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_GaSb_PmaxTMscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, GaSb_Ti_PmaxTMscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_GaSb_aoboptTEabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_GaSb_boloptTEabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, GaSb_Ti_aoboptTEabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, GaSb_Ti_boloptTEabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('Ti/GaSb max P^{TE}_{abs}/P^{TE}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_GaSb_aoboptTMabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_GaSb_boloptTMabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, GaSb_Ti_aoboptTMabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, GaSb_Ti_boloptTMabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_GaSb_aoboptTEscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_GaSb_boloptTEscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, GaSb_Ti_aoboptTEscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, GaSb_Ti_boloptTEscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_GaSb_aoboptTMscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_GaSb_boloptTMscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, GaSb_Ti_aoboptTMscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, GaSb_Ti_boloptTMscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);

%---------------------------------------------------------------TITANIUM/GE
load('Ti_Ge_Data.mat');
load('Ge_Ti_Data.mat');
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_Ge_PmaxTEabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, Ge_Ti_PmaxTEabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('Ti/Ge max P^{TE}_{abs}/P^{TE}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_Ge_PmaxTMabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, Ge_Ti_PmaxTMabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_Ge_PmaxTEscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, Ge_Ti_PmaxTEscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_Ge_PmaxTMscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, Ge_Ti_PmaxTMscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_Ge_aoboptTEabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_Ge_boloptTEabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, Ge_Ti_aoboptTEabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, Ge_Ti_boloptTEabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('Ti/Ge max P^{TE}_{abs}/P^{TE}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_Ge_aoboptTMabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_Ge_boloptTMabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, Ge_Ti_aoboptTMabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, Ge_Ti_boloptTMabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_Ge_aoboptTEscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_Ge_boloptTEscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, Ge_Ti_aoboptTEscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, Ge_Ti_boloptTEscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_Ge_aoboptTMscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_Ge_boloptTMscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, Ge_Ti_aoboptTMscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, Ge_Ti_boloptTMscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);

%-------------------------------------------------------------TITANIUM/INAS
load('Ti_InAs_Data.mat');
load('InAs_Ti_Data.mat');
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_InAs_PmaxTEabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, InAs_Ti_PmaxTEabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('Ti/InAs max P^{TE}_{abs}/P^{TE}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_InAs_PmaxTMabs, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, InAs_Ti_PmaxTMabs, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_InAs_PmaxTEscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, InAs_Ti_PmaxTEscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_InAs_PmaxTMscat, '-bo', 'LineWidth', 2.1);
plot(lambda01/nm, InAs_Ti_PmaxTMscat, '-mo', 'LineWidth', 1.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
xlim([400, 700]); ylim([1, 3]);
figure;
subplot(2, 2, 1);
hold;
plot(lambda01/nm, Ti_InAs_aoboptTEabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_InAs_boloptTEabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, InAs_Ti_aoboptTEabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, InAs_Ti_boloptTEabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('Ti/InAs max P^{TE}_{abs}/P^{TE}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 2);
hold;
plot(lambda01/nm, Ti_InAs_aoboptTMabs, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_InAs_boloptTMabs, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, InAs_Ti_aoboptTMabs, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, InAs_Ti_boloptTMabs, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{abs}/P^{TM}_{0,abs}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 3);
hold;
plot(lambda01/nm, Ti_InAs_aoboptTEscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_InAs_boloptTEscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, InAs_Ti_aoboptTEscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, InAs_Ti_boloptTEscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TE}_{scat}/P^{TE}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);
subplot(2, 2, 4);
hold;
plot(lambda01/nm, Ti_InAs_aoboptTMscat, '-go', 'LineWidth', 2.1);
plot(lambda01/nm, Ti_InAs_boloptTMscat, '-ro', 'LineWidth', 2.1);
plot(lambda01/nm, InAs_Ti_aoboptTMscat, '-co', 'LineWidth', 2.1);
plot(lambda01/nm, InAs_Ti_boloptTMscat, '-mo', 'LineWidth', 2.1);
xlabel('\lambda_0 (nm)');
title('max P^{TM}_{scat}/P^{TM}_{0,scat}');
legend('a/b', 'b/\lambda_0');
grid;
xlim([lambda0min/nm, lambda0max/nm]);