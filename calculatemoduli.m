clc;
% Process BM data
pre = input('Delay1 (s): ');
post = input('Delay2 (s): ');
pulse = input('Pulse (s): ');
analyzedind = input('Indices to analyze: ');
analyzeyn = 1;

% Bimorph properties
KD = 2.48e-3; % 1/(V?m)
KF = 2.62e-3; % N/V
CB = 1.551e-3; % m/N
KB = 1/CB;    % N/m
FB0 = 7.74e-4;% N/V
DU0 = 1.20032e-6;% m/V
MB = 0.386e-3; % kg
fr = 1 / (2*pi*sqrt(MB*CB));    % Hz
CP = 8.9e-9;      % F
tau = 1/fr; % s
    % Displacement of bimorph is Db = Cb*Fb = Fb/Kb
    % Force on a bimorph with standard size -> Fb0 = 6.81e-4 N/V, 
    % or Fb = Eb * Kf * w / l
    % With an elastic load, D1 = Db * (C1 / (C1 + Cb))
    % and F1 = Fb * (Cb / (Cb + C1))
    
clear Xd_split
for ind = analyzedind
    disp(['Analyzing index ' num2str(ind)]);
    %{
    Nf(ind) = input('Number of pulses: ');
    fmin(ind) = input('Minimum position (µm): ');
    stepsize(ind) = input('Step size (µm): ');
    %}
    
    Nf(ind) = logdata.data(ind,14);
    fmin(ind) = logdata.data(ind,12);
    stepsize(ind) = logdata.data(ind,13);
    basedisplacements{ind} = fmin(ind):stepsize(ind):(stepsize(ind)*(Nf(ind)-1)+fmin(ind)); 
    
    pulselength = length(Xd{ind})/Nf(ind);
    if mod(pulselength,1) == 0
    for j = 1:Nf(ind)
        Xd_split{ind}(:,j) = Xd{ind}(1+(j-1)*pulselength:j*pulselength);
        meanpre(j,ind) = mean(Xd_split{ind}(pre*Fs-round(pre*Fs/20):pre*Fs+round(pre*Fs/20),j));
        pulsemax(j,ind) = max((Xd_split{ind}(pre*Fs:(pre+pulse)*Fs,j)));
        pulsemin(j,ind) = min((Xd_split{ind}(pre*Fs:(pre+pulse)*Fs,j)));
        qr = find(Xd_split{ind}(:,j) == pulsemax(j));
        if isempty(qr) == 0
            maxind(j,ind) = qr(1);
        else
            maxind(j,ind) = 6;
        end
        pulsemean(j,ind) = mean(Xd_split{ind}((maxind(j,ind)-3:maxind(j,ind)+3),j));
        Xd_split_integral{ind}(:,j) = cumtrapz(tvec{ind}(1:length(Xd_split{ind}(:,j))),Xd_split{ind}(:,j));
        pulse_int_max(j,ind) = max(Xd_split_integral{ind}((pre+pulse)*Fs:1.2*(pre+pulse)*Fs,j));
    end
    end
end

%% Plot
close all
clc
area0(1:5)=10;
% length0(12:14) = length0(6:8);
tstart = 0; % 1 + tstart
tend = 0;   % end - tend

controlmoduli = [1:5];
expmoduli = [6:9];
setfiguredefaults(length(controlmoduli) + length(expmoduli))

for ind = analyzedind
    forceapplied{ind} = pulsemax(1:length(basedisplacements{ind}),ind)'.*FB0;
    deflections{ind} = pulsemax(1:length(basedisplacements{ind}),ind)'.*DU0;
%     indentations{ind} = displacements{ind} - forceapplied{ind}./KB;
    indentations{ind} = basedisplacements{ind}.*1e-6 - basedisplacements{ind}(1).*1e-6 - deflections{ind};
    stresses{ind} = forceapplied{ind} ./ area0(ind);
%     strains{ind} = ( (basedisplacements{ind}.*1e-6-basedisplacements{ind}(1)*1e-6) - displacements{ind}) ./ length0(ind);
    strains{ind} = indentations{ind} ./ length0(ind) ;
    figure(1);subplot(1,2,1);
    plot(strains{ind},stresses{ind});hold all;grid on;
    xlabel('Strain (dL/L)');ylabel('Stress (F/A) (N?m^-^2)')
    figure(1);subplot(1,2,2);
    plot(indentations{ind}.*1e6,forceapplied{ind});hold all;grid on;
    xlabel('Indentation (µm)');ylabel('Force (N)')
    [fits{ind} gofs{ind}] = fit(strains{ind}(1 + tstart:end - tend)',stresses{ind}(1 + tstart:end - tend)','poly1');
    [fits2{ind} gofs2{ind}] = fit(indentations{ind}(1 + tstart:end - tend)',forceapplied{ind}(1 + tstart:end - tend)','poly1');
    modulus(ind) = fits{ind}.p1;
    stiffnesses(ind) = fits2{ind}.p1;
    R2(ind) = gofs{ind}.rsquare;
    legends{ind} = comments{ind}(1);
%     plot(displacements{ind},pulsemax(1:length(basedisplacements{ind}),ind).*FB0);hold all;
end


meanmodulus_ctrl = mean(modulus(controlmoduli));
semmodulus_ctrl = std(modulus(controlmoduli))./sqrt(length(controlmoduli));
meanmodulus_exp = mean(modulus(expmoduli));
semmodulus_exp = std(modulus(expmoduli))./sqrt(length(expmoduli));

% Calculate elastic modulus of material, assuming it is in series with all
% other elements
material_modulus = 1/( 1/meanmodulus_exp - 1/meanmodulus_ctrl );    % series
material_modulus = meanmodulus_ctrl - meanmodulus_exp;    % parallel

material_modulus_sem = material_modulus * sqrt((semmodulus_ctrl./meanmodulus_ctrl)^2 + (semmodulus_exp./meanmodulus_exp)^2);
figure(1);subplot(1,2,1);
disp(['Material modulus: ' num2str(material_modulus) ' ± ' num2str(material_modulus_sem) ' Pa'])
title(['Material modulus: ' num2str(material_modulus) ' ± ' num2str(material_modulus_sem) ' Pa'])
legend(legends);

%% Calculate various moduli
inds{1} = [1:3];
inds{2} = [4:6];
inds{3} = [7:9];
inds{4} = 12:14;
inds{5} = [15 18];
inds{6} = [19 20 23];
inds{7} = 25:27;
inds{8} = 28:29;
inds{9} = 32:34;
inds{10} = 36:39;
inds{11} = 40:43;
inds{12} = 45:49;


area0 = [
7.03E-05
5.67E-05
7.13E-05
5.63E-05
7.67E-05
5.57E-05
8.27E-05
7.87E-05
7.33E-05
7.90E-05
7.57E-05
8.47E-05
];
length0 = [
2.84E-03
3.53E-03
2.80E-03
3.55E-03
2.61E-03
3.59E-03
2.42E-03
2.54E-03
2.73E-03
2.53E-03
2.64E-03
2.36E-03
];

indexkey{1} = '2% - 1';
indexkey{2} = '4% - 1';
indexkey{3} = '2% - 2';
indexkey{4} = '4% - 2';
indexkey{5} = '2% - 3';
indexkey{6} = '4% - 3';
indexkey{7} = '1% - 1';
indexkey{8} = '0.5% - 1';
indexkey{9} = '1% - 2';
indexkey{10} = '0.5% - 2';
indexkey{11} = '1% - 3';
indexkey{12} = '0.5% - 3';
colors = {'r--.','g--.','r--o','g--o','r--*','g--*','c--.','b--.','c--o','b--o','c--*','b--*'};
tstart = 0; % 1 + tstart
tend = 0;   % end - tend
for m = 1:length(area0)
for ind = inds{m}
forceapplied{ind} = pulsemax(1:length(basedisplacements{ind}),ind)'.*FB0;
deflections{ind} = pulsemax(1:length(basedisplacements{ind}),ind)'.*DU0;
indentations{ind} = basedisplacements{ind}.*1e-6 - basedisplacements{ind}(1).*1e-6 - deflections{ind};
stresses{ind} = forceapplied{ind} ./ area0(m);
strains{ind} = indentations{ind} ./ length0(m) ;
figure(1);subplot(1,2,1);
plot(strains{ind},stresses{ind},colors{m});hold all;grid on;
xlabel('Strain (dL/L)');ylabel('Stress (F/A) (N?m^-^2)')
figure(1);subplot(1,2,2);
plot(indentations{ind}.*1e6,forceapplied{ind},colors{m});hold all;grid on;
xlabel('Indentation (µm)');ylabel('Force (N)')
[fits{m}{ind} gofs{m}{ind}] = fit(strains{ind}(1 + tstart:end - tend)',stresses{ind}(1 + tstart:end - tend)','poly1');
%[fits2{m}{ind} gofs2{m}{ind}] = fit(indentations{ind}(1 + tstart:end - tend)',forceapplied{ind}(1 + tstart:end - tend)','poly1');
modulus{m}(ind) = fits{m}{ind}.p1;modulus{m}(modulus{m}==0)=NaN;
%stiffnesses{m}(ind) = fits2{m}{ind}.p1;
R2{m}(ind) = gofs{m}{ind}.rsquare;R2{m}(R2{m}==0)=NaN;
end
meanmodulus(m) = nanmean(modulus{m});
meanR2(m) = nanmean(R2{m});
end
%% Manual curve fitting
colors = {'r--.','g--.','r--o','g--o','r--*','g--*','c--.','b--.','c--o','b--o','c--*','b--*'};
for m = 1:length(area0)
    for ind = inds{m}
        disp(['m = ' num2str(m) '; ind = ' num2str(ind)])
        yintercept{m}(ind) = fits{m}{ind}.p2;
        if length(yintercept{m}) >= ind
        if isempty(yintercept{m}(ind))
            yintercept{m}(ind) = fits{m}{ind}.p2;
        end
        else
            yintercept{m}(ind) = fits{m}{ind}.p2;
        end
        if R2{m}(ind) < 0.95
            cftool(strains{ind},stresses{ind});
            modulus{m}(ind) = input('Slope (p1): ');modulus{m}(modulus{m}==0)=NaN;
            yintercept{m}(ind) = input('Y intercept (p2): ');
            R2{m}(ind) = input('R2: ');R2{m}(R2{m}==0)=NaN;
        end
        
        xintercept{m}(ind) = -yintercept{m}(ind) / modulus{m}(ind);
        figure(2);
        plot(strains{ind}-xintercept{m}(ind),stresses{ind},colors{m});
        exclxstart = input('Exclude how many from start: ');
        exclxend = input('Exclude how many from end: ');
        close(2);
        figure(1);subplot(1,2,1);
        plot(strains{ind}(1+exclxstart:end-exclxend)-xintercept{m}(ind),stresses{ind}(1+exclxstart:end-exclxend),colors{m});hold all;grid on;
        xlabel('Strain (dL/L)');ylabel('Stress (F/A) (N?m^-^2)')
        figure(1);subplot(1,2,2);
        plot(indentations{ind}(1+exclxstart:end-exclxend).*1e6,forceapplied{ind}(1+exclxstart:end-exclxend),colors{m});hold all;grid on;
        xlabel('Indentation (µm)');ylabel('Force (N)')
    end
    
    meanmodulus(m) = nanmean(modulus{m});
    semmodulus(m) = nanstd(modulus{m}) / sqrt(length(inds{m}));
    meanR2(m) = nanmean(R2{m});
end

%%
ctrl_ind = 1;
exp_ind = 2:5;
q = 1;
for m = 1:length(area0)
    material_modulus(m) = 1/( 1/meanmodulus(m) - 1/(meanmodulus(1)) ); 
    material_modulus_sem(m) = material_modulus(m) * sqrt((semmodulus(1)./meanmodulus(1))^2 + (semmodulus(m)./meanmodulus(m))^2);
end
