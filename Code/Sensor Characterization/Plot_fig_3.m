%% Code by Natalia Sanchez-Tamayo
% Date: March 5th, 2024 
% Last updated: March 10th, 2026
% This script plots Figure 3, which characterizes the sensor array
%% Import data of intentations. On taxel (T3) indentations
force_array_T3=readmatrix("Cyclic measurement of sensor_Data_25_1_24_3.txt");
cap_raw_T3=readmatrix("Sensor_response_cyclic_3.txt");

%% 3.A1 
T3_idx=7;
% Linear calibration coefficients (ADC counts → pF)
% Derived from least-squares fit to calibration data
b=16.904751517474693; 
m=-0.0017790184866274645;

%average baseline across taxels is 
baseline_c=b+m*9362;
time_step=0.02;
% process and sync data %I used med filter by 5 to plot the indentations
[force_T3,cap_T3,cycle_T3,time_T3]=sync_data(force_array_T3,cap_raw_T3,T3_idx,1.16,"Taxel 3",m,b,time_step)

%% Function: sync_data
% This function synchronizes all capacitance values and cycle numbers to the force measurements.
function [force_f,cap_f,cycle_f,time_sync]=sync_data(force,cap,push_id,time_offset,name,m,b,time_step)
    
    %cap delay
    mt=0.0005
    bt=0.5016

    % this function syncs all capacitance valuesl and cycle number to the forces
    xq = 0:pi/16:2*pi;
    %original forces are in 100 ms, 0.1 and cap close to 60hz 16ms, 0.015
    step=time_step
    time_sync=0:step:force(end,1); % finer grid of time to interpolate to 
    cap_time=cap(:,end)/1000-time_offset; 
    cap_f=zeros(length(time_sync),16);


    
    cap_time=cap_time+cap_time*mt; %adjusting for delay using only slope

    figure(90)
    yyaxis left
    plot(cap_time(:),cap(:,push_id))
    hold on
    yyaxis right
    hold on
    plot(force(:,1),force(:,2))


    % filter all capacitance values by median filter
    for j=1:16
        if j==7
            fff=6;
        end
        cap(:,j)=medfilt1(cap(:,j),3); 
        cap_f(:,j) = -interp1(cap_time,cap(:,j),time_sync)*(-m); 
    end 
    
    %     interpolate values from force and cycle 
        force_f = interp1(force(:,1), force(:,2),time_sync);
        cycle_f=interp1(force(:,1),force(:,6),time_sync);
        

        close all
        fig = figure(5);
        left_color = [224/255 119/255 119/255];
        right_color = [18/255 17/255 13/255];
        ored=[231/255 111/255 81/255];
        ora=[244/255 162/255 97/255];
        hold on 
        yyaxis right
        ylabel('Capacitance [pF]')
        plot(time_sync,cap_f(:,push_id),'color',left_color);
        set(gca,'YColor',left_color);
        hold on
        plot(time_sync,cap_f(:,push_id-1),'-','color',ored);
        plot(time_sync,cap_f(:,push_id+1),'-','color',ora);
        ylim([-0.18 2.25])

        yyaxis left
        ylabel('Force [N]')
        plot(time_sync,force_f,'color',right_color)
        set(gca,'YColor',right_color);
        ylim([-0.25 3.25])

        % title(name)
        x0=100;
        y0=100;
        width=850;
        height=160;
        set(gcf,'position',[x0,y0,width,height])

        % CLOSE UP
        fig = figure(50);
        left_color = [73/255 171/255 137/255];
        right_color = [137/255 6/255 32/255];
        hold on 
        yyaxis right
        ylabel('Capacitance [pF]')
        plot(time_sync,cap_f(:,push_id),'color',left_color);
        set(gca,'YColor',left_color);
        ylim([-200*(-m) 1100*(-m)])
        hold on
        yyaxis left
        ylabel('Force [N]')
        plot(time_sync,force_f,'color',right_color)
        set(gca,'YColor',right_color);
        ylim([-0.25 2.1])
        xlim([300,570]);
        % title(name)
        x0=100;
        y0=100;
        width=350;
        height=250;
        set(gcf,'position',[x0,y0,width,height])

        figure(21)
        TF = islocalmin(force_f,'MinProminence',0.04);
        Tcap = islocalmin(cap_f(:,push_id),'MinProminence',0.5);
        plot(time_sync,force_f,time_sync(TF),force_f(TF),'r*')
        figure(22)
        plot(time_sync,cap_f(:,push_id),time_sync(Tcap),cap_f(Tcap,push_id),'r*')
end

%%

left_color = [224/255 119/255 119/255];
right_color = [18/255 17/255 13/255];
name_vid="Indent_T3_v5";
sampleRate=1/(time_T3(2)-time_T3(1));

figure('Position',[100 100 800 600]); 
count=1;

cap_T3(:,:)=cap_T3(:,:)-mean(cap_T3(1:100,1));
for i = 1:length(time_T3)/8
    count=count+1;
    disp(i)
    subplot(2,1,1);
    cla; 
    h = plot(time_T3(1:i), cap_T3(1:i,:)*1000,'LineStyle',':','LineWidth',1);
    h(T3_idx).Color = left_color;
    h(T3_idx).LineWidth = 2;
    h(T3_idx).LineStyle='-';
    uistack(h(T3_idx), 'top');

    ylabel('C - C_{0} [fF]', 'FontSize', 15, 'FontWeight', 'bold');
    xlabel('Time [s]', 'FontSize',15, 'FontWeight', 'bold');
    xlim([0, 175]);
    ax = gca; 
    ax.FontSize = 15; 

    subplot(2,1,2);
    plot(time_T3(1:i), force_T3(1:i)*1000,'color',right_color,'LineWidth',2);
    ylabel('Force [mN]', 'FontSihze', 15, 'FontWeight', 'bold');
    xlabel('Time [s]', 'FontSize', 15, 'FontWeight', 'bold');
    ax = gca; 
    ax.FontSize = 15; 
    xlim([0, 175]);

    drawnow
    F(i) = getframe(gcf);
end

v = VideoWriter(name_vid);
v.FrameRate = sampleRate;
open(v);
writeVideo(v, F);
close(v);

%% 3.B — Minimum Detectable Force
% Plots the sensor response at low force levels to characterize
% the noise floor and detection threshold of Taxel 3.
% fig = figure;
left_color = [73/255 171/255 137/255];
right_color = [4/255 42/255 98/255];

figure (100)
step=0.1; 
yyaxis left
plot(time_T3(:),force_T3(:),'LineWidth',0.8,'color',right_color);
set(gca,'YColor',right_color);

ylim([-0.009 0.020])
ylabel('Force (N)')
xlabel('Time [s]');
yyaxis right
ored=[231/255 111/255 81/255];
ora=[244/255 162/255 97/255];
plot(time_T3(:),cap_T3(:,T3_idx-1),'color',ored);
hold on
plot(time_T3(:),cap_T3(:,T3_idx+1),'color',ora);
hold on
plot(time_T3(:),cap_T3(:,T3_idx),'color',left_color);

xlim([0 50])
ylim([-0.004 0.008])

set(gca,'YColor',left_color);
x0=100;
y0=100;
width=350;
height=250;
set(gcf,'position',[x0,y0,width,height])

%% min detectable force alignes
L1=0; L2=10; R1=0; R2=0.0035; pad=1;
a=(R2-R1)/(L2-L1); b=R1-a*L1;
Lmin=min(L1,L2); Lmax=max(L1,L2); span=Lmax-Lmin;
yl_left=[Lmin-pad*span, Lmax+pad*span];

%Plot the minimum detectable force 
right_color =[0/255 0/255 0/255];
left_color =  [159/255 18/255 40/255];

figure (100)
step=0.1;
yyaxis left
plot(time_T3(:),force_T3(:)*1000,'LineWidth',0.8,'color',right_color);
set(gca,'YColor',right_color);

ylim([-0.009 0.020])
ylabel('Force [mN]')
xlabel('Time [s]');
yyaxis right

blue1 = [96/255 202/255 204/255];
blue2 = [80/255 128/255 186/255];

plot(time_T3(:),cap_T3(:,T3_idx-1),'color',blue1);
hold on
plot(time_T3(:),cap_T3(:,T3_idx+1),'color',blue2);
hold on
plot(time_T3(:),cap_T3(:,T3_idx),'color',left_color);

xlim([0 50])

ylabel('Capacitance [pF]')
set(gca,'YColor',left_color);
x0=100;
y0=100;
width=350;
height=250;
set(gcf,'position',[x0,y0,width,height])


yyaxis left,  ylim(yl_left)
yyaxis right, ylim(a*yl_left + b) 

%%
%Plot closeup 2
right_color =[181/255 146/255 160/255];
left_color =  [103/255 0/255 31/255];
blue1 = [164/255 206/255 227/255];
blue2 = [34/255 116/255 165/255];


figure
step=0.1; 
subplot(2,1,1)
plot(time_T3(:),force_T3(:),'LineWidth',0.8,'color',right_color);
set(gca,'YColor',right_color);
ylim([0.005 0.05])
ylabel('Force (N)')
xlabel('Time [s]');
xlim([800 990])
% yyaxis left
subplot(2,1,2)
plot(time_T3(:),cap_T3(:,T3_idx),'color',left_color);
xlim([800 990])
ylim([0.01 0.05])
ylabel('Capacitance [pF]')

set(gca,'YColor',left_color);
x0=100;
y0=100;
width=350;
height=250;
set(gcf,'position',[x0,y0,width,height])



%superimposed
figure (10)
step=0.1;
yyaxis right
plot(time_T3(:),force_T3(:).*1000,'LineWidth',0.8,'color',right_color);
hold on
yline(10)
yline(0)
set(gca,'YColor',right_color);
ylim([-9 20])
ylabel('Force [mN]')
xlabel('Time [s]');

xlim([0 400])
yyaxis left


plot(time_T3(:),cap_T3(:,T3_idx).*1000,'color',left_color);

xlim([0 50])
ylim([-4 8.9])
ylabel('Capacitance [fF]')
set(gca,'YColor',left_color);
x0=100;
y0=100;
width=350;
height=250;
set(gcf,'position',[x0,y0,width,height])

figure (11)
step=0.1;
yyaxis left
plot(time_T3(:),force_T3(:).*1000,'LineWidth',0.8,'color',right_color);
set(gca,'YColor',right_color);
ylim([-10 40])
ylabel('Force [mN]')
xlabel('Time [s]');
xlim([0 400])
yline(10)
yline(0)
yyaxis right


plot(time_T3(:),cap_T3(:,T3_idx-1).*1000,'color',blue1);
hold on
plot(time_T3(:),cap_T3(:,T3_idx+1).*1000,'color',blue2);
hold on
plot(time_T3(:),cap_T3(:,T3_idx).*1000,'color',left_color);
xlim([800 980])
ylim([-20 60])
ylabel('Capacitance [fF]')
set(gca,'YColor',left_color);
x0=100;
y0=100;
width=350;
height=250;
set(gcf,'position',[x0,y0,width,height])

L1=0; L2=10.24; R1=0; R2=14.23; pad=3.5;
a=(R2-R1)/(L2-L1); b=R1-a*L1;
Lmin=min(L1,L2); Lmax=max(L1,L2); span=Lmax-Lmin;
yl_left=[Lmin-pad/2*span, Lmax+pad*span];

yyaxis left,  ylim(yl_left)
yyaxis right, ylim(a*yl_left + b)


%%
%Peaks
left_color = [73/255 171/255 137/255];
right_color = [4/255 42/255 98/255];

%superimposed
figure
step=0.1; 
yyaxis right
plot(time_T3(:),force_T3(:).*1000,'LineWidth',0.8,'color',right_color);
hold on
yline(10)
yline(0)
set(gca,'YColor',right_color);
ylim([-9 20])
ylabel('Force [mN]')
xlabel('Time [s]');

xlim([0 400])
yyaxis left

% subplot(2,1,2)
plot(time_T3(:),cap_T3(:,T3_idx).*1000,'color',left_color);
xlim([0 50])
ylim([-4 8.9])
ylabel('Capacitance [fF]')
set(gca,'YColor',left_color);
x0=100;
y0=100;
width=350;
height=250;
set(gcf,'position',[x0,y0,width,height])

figure
step=0.1;
yyaxis right
plot(time_T3(:),force_T3(:).*1000,'LineWidth',0.8,'color',right_color);
set(gca,'YColor',right_color);
ylim([2970 3010])
ylabel('Force [mN]')
xlabel('Time [s]');
xlim([0 400])
yline(10)
yline(0)
yyaxis left

plot(time_T3(:),cap_T3(:,T3_idx).*1000,'color',left_color);
xlim([890 1100])
ylim([1790 1800])
ylabel('Capacitance [fF]')
set(gca,'YColor',left_color);
x0=100;
y0=100;
width=350;
height=250;
set(gcf,'position',[x0,y0,width,height])


%% 3.C plot SNR
figure
plot_noise(force_T3,cap_T3,cycle_T3,time_T3,T3_idx)
% SNR is plotted with ulfiltered data for the first indentation during
% loading. ploted with 60 samples per point.

function[]=plot_noise(force_f,cap_f,cycle_f,time_sync,t_id)

    left_color = "#175290";
    right_color = [4/255 42/255 98/255];

    num_c=60;
    cap_noise=cap_f(1:num_c,t_id);
    
    % Plot SNR for the first indentations:
    idx_end=5218; % First peak
    list_idx=1:num_c:idx_end;
    tras=zeros(1,length(list_idx));
    forces=zeros(1,length(list_idx));

    for i=list_idx
        signal=cap_f(i:i+num_c-1,t_id);
        tras(i)=snr(signal,cap_noise);
        forces(i)=mean(force_f(i:i+num_c-1));
    end
    scatter(forces,tras,90,'.','MarkerEdgeColor',left_color,'MarkerFaceColor',left_color)
    % set(gca,'YColor',left_color);
    x0=100;
    y0=100;
    xlim([0.01,3])
    ylim([0, 70])
    width=350;
    height=250;
    set(gcf,'position',[x0,y0,width,height])
    set(gca, 'XScale', 'log')
    ylabel('Signal to noise ratio [dB]')
    xlabel('Force [N]')

end

%% 
% 3.D Plot the sensor response for 10 cycles 
plot_repeatibility(force_T3,cap_T3,cycle_T3,time_T3,T3_idx)
% cap_noise=cap_T3H(1:200,T3_idx);

function []=plot_repeatibility(force_f,cap_f,cycle_f,time_sync,T3_idx)

    % For each cycle find the values of capacitance during loading that
    % reach a specific value, for x number of force:
    % Index of the peak for the first full loading cycle (manually verified)
    peaks_idx=[5218,11835,18425,24978,31525,38056,44575,51083,57556,64030];

    len_p=3200;
    forces_bins=0.5:0.5:3.5;
    capacitances=zeros(length(peaks_idx),length(forces_bins));
    noise=zeros(length(forces_bins));

    for i=1:length(forces_bins)
        for j=1:length(peaks_idx)
               x = force_f(peaks_idx(j):peaks_idx(j)+len_p); % extract two arryays for each indentation
               y = cap_f(peaks_idx(j):peaks_idx(j)+len_p,T3_idx)' ;
               name=sprintf('Cycle %.1f',j)
               hold on
               iwant = y(x>forces_bins(i)-0.004 & x <forces_bins(i)+0.004) ; % all cap values when f is within 0.003 of the desired value
               capacitances(j,i)=mean(iwant);
               if j==1
                   cap_signal=y;
               end
        end

    end
   

    figure(40)
    ylabel('\drift sensor response per cycles')
    xlabel('Force [N]')

    cycles=1:1:10;

    for k=1:length(forces_bins)-1
        disp(k)
        hold on
        name=sprintf('%.1f N',forces_bins(k))
       plot(cycles(:),capacitances(:,k),'-o','MarkerSize',4, 'DisplayName',name);
    end

    % Mean across repetitions (rows)
    mean_per_force = mean(capacitances*1000, 1);
    
    % Standard deviation across repetitions (rows)
    std_per_force = std(capacitances*1000, 0, 1);   % 0 -> normalize by (n-1) for sample std

    legend('Orientation','vertical','NumColumns', 1)
    ylim([0 2])
    xlabel('Cycle number')
    ylabel('C - C_{0} [pF]')

    x0=100;
    y0=100;
    width=350;
    height=250;
    set(gcf,'position',[x0,y0,width,height])
    figure (60)
    force_levels = 1:size(capacitances,2);  
    errorbar(force_levels, mean_per_force, std_per_force, 'o-');
    xlabel('Force Level');
    ylabel('Capacitance');
    title('Capacitance vs Force Level (Mean ± Std)');
    grid on;

    figure(61);
    errorbar(force_levels, mean_per_force, std_per_force, 'o-','LineWidth',1.5);
    xlabel('Force Level');
    ylabel('Capacitance');
    title('Capacitance vs Force Level (Zoomed Inset)');
    grid on;
    
    % Add inset to zoom in on error bars
    axes('Position',[.6 .6 .25 .25]);
    box on;
    errorbar(force_levels, mean_per_force, std_per_force, 'o-','LineWidth',1.5);
    xlim([min(force_levels)-0.5, max(force_levels)+0.5]);
    ylim([min(mean_per_force)-2*max(std_per_force), max(mean_per_force)+2*max(std_per_force)]);

    cv_per_force = (std_per_force ./ mean_per_force) * 100;

    figure(63);
    bar(forces_bins(1:end-2), cv_per_force(1:end-2), 'FaceColor',[0.8 0.4 0.2]);
    xlabel('Force [mN]');
    ylabel('Coefficient of Variation (%)');
    title('Repeatibility');
    legend('Relative Variability of Capacitance');
    xlim([0, 3.5]);
end


%%
close all
plot_hysteresis(force_T3,cap_T3,cycle_T3,time_T3,T3_idx)
% cap_noise=cap_T3H(1:200,T3_idx);


function []=plot_hysteresis(force_f,cap_f,cycle_f,time_sync,T3_idx)

    peaks_idx=[5218,11835,57556,64030];

    % %only the cycle 1 and 10. full cycles
    peaks_idx=[5218,64030];

    len_p=3200;%3301;
    figure(20)
    forces_bins=0.5:0.5:3.5;
    capacitances=zeros(length(peaks_idx),length(forces_bins));
    noise=zeros(length(forces_bins));

    for j=1:length(peaks_idx)
           x = force_f(peaks_idx(j)-len_p:peaks_idx(j)+len_p); % extract two arryays for each indentation
           y = cap_f(peaks_idx(j)-len_p:peaks_idx(j)+len_p,T3_idx)' ;
           tt = time_sync(peaks_idx(j)-len_p:peaks_idx(j)+len_p); 
           % y_0 = y(200-end:end) ;
           name=sprintf('Cycle %.1f',j)
           hold on
           figure(20)
           plot(x,medfilt1(y,3),'DisplayName',name)
           figure(21)
           plot(x,medfilt1(y,3),['DisplayName'],name)
    end
    figure(20)
    x0=100;
    y0=100;
    width=350;
    height=250;
    set(gcf,'position',[x0,y0,width,height])
    legend('cycle 1','cycle 2',' cycle 9','cycle 10')

    ylabel('C - C_{0} [pF]')
    xlabel('Time [s]')

    figure(21)
    xlim([0.5,0.6])
    x0=100;
    y0=100;
    width=100;
    

    height=100;
    set(gcf,'position',[x0,y0,width,height])

   
end

%% Sensitivity

close all
[force_T3s,cap_T3s,cycle_T3s,time_T3s]=sync_data(force_array_T3,cap_raw_T3,T3_idx,1.16,"Taxel 3",m,b,0.02)
%
[sensitivity]=plot_sensitivity(force_T3s,cap_T3s,cycle_T3s,time_T3s,T3_idx,m);


function [sensitivity]=plot_sensitivity(force_f,cap_f,cycle_f,time_sync,T3_idx,m)

    peaks_idx=[5218];

    len_p=4100;

    j=1;
    x = force_f(peaks_idx(j)-len_p:peaks_idx(j)); % extract two arryays for each indentation
    y = cap_f(peaks_idx(j)-len_p:peaks_idx(j),T3_idx)' ;
    tt = time_sync(peaks_idx(j)-len_p:peaks_idx(j)); 

   
    F = sgolayfilt(x,3,15);
    C = sgolayfilt(y,3,15);
   

    delta=0.01;% for specific delta force 20 gr
    FQ = 0.0:delta:3; % Points to query the force

    CQ=interp1(F,C,FQ); % queried Capacitnace values.

    % calculate the delta forces and c at each step
    delta_f=FQ(2:end)-FQ(1:end-1);
    delta_c=CQ(2:end)-CQ(1:end-1);
    % 
    sensitivity=delta_c./delta_f;

    % scatter(FQ(2:end),sensitivity)
    ylim([0,1.2])

    xc=FQ(2:end);
    yc=sensitivity;

    fitted_model = fit(xc', yc','poly9');
    n=feval(fitted_model, xc);
    residuals = yc - n';

    bins=delta*10;
    bin_f = 0.01:bins:3+bins*1;

    f_res = zeros(length(bin_f)-1, 2);
    s_res = zeros(length(bin_f)-1, 2);

    for j = 2:length(bin_f)
        res_list = find(xc < bin_f(j) & xc > bin_f(j-1));
        f_res(j,1) = mean(xc(res_list));
        s_res(j,2) = std(yc(res_list)'-feval(fitted_model, xc(res_list))); % std of the residuals for all xc points between these force ranges
        f_res(j,2) = std(xc(res_list));
        s_res(j,1) = mean(feval(fitted_model, xc(res_list)));
    end

    % Plot sensitivity
    figure (12);
    std_res = sgolayfilt(s_res(:,2),3,9);
    x_plot=f_res(:,1);
    y_plot=s_res(:,1);

    fill_color= [82, 119, 122]/256;
    fill([x_plot; flipud(x_plot)], [y_plot-std_res; flipud(y_plot+std_res)], fill_color, 'FaceAlpha', 0.2);
    hold on
    plot(xc, feval(fitted_model, xc), '-k')
    % xlabel('Force [N]', 'FontSize', 15);
    % ylabel('Sensitivity [pF/N]', 'FontSize', 15);
    set(gca, 'FontSize', 8);
    xlim([0.01, 3]);
    ylim([-0.1, 1.2]);
    x0=100;
    y0=100;
    width=350;
    height=250;
    set(gcf,'position',[x0,y0,width,height])

end





