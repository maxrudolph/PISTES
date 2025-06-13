
clear;
close all;
% run the thermal evolution model
seconds_in_year = 3.1558e7;
parameters.no_stress_time = 3.5e9;
% parameters.crust_heat_fraction=0.5;
% parameters.viscosity = 3e20;
parameters.do_plots = false;

nvisc=10;
visc = logspace(20,22,nvisc);
nhf=11;
hf = linspace(0.3,0.7,nhf);
nstress = 10;
nst = linspace(0,3.5e9,nstress);

allresults = cell(nvisc,nhf);
for ivisc=1:nvisc
    parfor ihf=1:nhf
        p = parameters;
        p.viscosity = visc(ivisc);
        p.crust_heat_fraction = hf(ihf);
        allresults{ivisc,ihf} = mars_thermal_evolution_and_stress(p);
    end
end

save(strcat('all_results_visc_crusthf', string(datetime), '.mat'),'-v7.3')

clear allresults;

for istress=1:nstress
    parfor ihf=1:nhf
        p = parameters;
        p.viscosity = 3e20;
        p.no_stress_time = nst(istress);
        p.crust_heat_fraction = hf(ihf);
        allresults{istress,ihf} = mars_thermal_evolution_and_stress(p);
    end
end

save(strcat('all_results_stress_crusthf', string(datetime), '.mat'),'-v7.3')


%% Make some plots to explore range from all models
% plot the crossover depth at the end of the calculation
max_stress_depth = cellfun( @(x) x.maximum_stress_depth(x.last_isave-1),allresults);
max_differential_stress = cellfun( @(x) x.maximum_differential_stress(x.last_isave-1),allresults);
min_differential_stress = cellfun( @(x) x.minimum_differential_stress(x.last_isave-1),allresults);
stress_crossover_depth = cellfun( @(x) x.stresss_crossover_depth(x.last_isave-1),allresults);
mantle_temperature=cellfun( @(x) x.Tm(x.last_isave-1),allresults);
final_lid_thickness=cellfun( @(x) x.state.Ro-x.state.Ri + x.z(x.last_isave-1),allresults);

f=figure()
f.Position(3:4) = [530   602];
t=tiledlayout(3,2,"TileSpacing","compact","Padding","none");
nexttile;
contourf(hf,visc,max_stress_depth/1e3);
set(gca,'YScale','log')
hcb = colorbar();
hcb.Label.String = 'Depth (km)';
set(gca,'PlotBoxAspectRatio',[1 1 1])

xlabel('Crustal heating fraction')
ylabel('\eta_0 (Pa-s)')
title("Depth of max. \sigma_t-\sigma_r")

nexttile
contourf(hf,visc,stress_crossover_depth/1e3);
set(gca,'YScale','log')
hcb=colorbar()
hcb.Label.String = 'Depth (km)';
xlabel('Crustal heating fraction')
ylabel('\eta_0 (Pa-s)')
title("Depth of \sigma_t-\sigma_r=0")
set(gca,'PlotBoxAspectRatio',[1 1 1])

nexttile
contourf(hf,visc,mantle_temperature);
set(gca,'YScale','log')
hcb=colorbar()
hcb.Label.String = 'T_m (K)';
xlabel('Crustal heating fraction')
ylabel('\eta_0 (Pa-s)')
title("Mantle temperature")
set(gca,'PlotBoxAspectRatio',[1 1 1])

nexttile
contourf(hf,visc,final_lid_thickness/1e3);
set(gca,'YScale','log')
hcb = colorbar();
hcb.Label.String = 'Thickness (km)';
xlabel('Crustal heating fraction')
ylabel('\eta_0 (Pa-s)')
title("Lid thickness")
set(gca,'PlotBoxAspectRatio',[1 1 1])

nexttile
contourf(hf,visc,max_differential_stress);
set(gca,'YScale','log')
hcb=colorbar();
hcb.Label.String = 'Stress (Pa)';
xlabel('Crustal heating fraction')
ylabel('\eta_0 (Pa-s)')
title("Max. \sigma_t-\sigma_r")
set(gca,'PlotBoxAspectRatio',[1 1 1])

nexttile
contourf(hf,visc,min_differential_stress);
set(gca,'YScale','log')
hcb=colorbar();
hcb.Label.String = 'Stress (Pa)';
xlabel('Crustal heating fraction')
ylabel('\eta_0 (Pa-s)')
title("Min. \sigma_t-\sigma_r")
set(gca,'PlotBoxAspectRatio',[1 1 1])
% exportgraphics(gcf,'tradeoff-eta-heating.pdf',"ContentType",'vector');

%% Plot the edge cases
for ivisc=[1 nvisc]
    for ihf=[1 nhf]
        results = allresults{ivisc,ihf};
        fieldnames = fields(results.state);
        for i=1:length(fieldnames)
            assignin("base",fieldnames{i},results.state.(fieldnames{i}));
        end

        %% Pseudocolor stress plot
        % mask = 1:(isave-1); % select only timesteps that exist
        mask = results.time <= results.last_save_time;
        isave = results.last_isave;
        save_depths = results.save_depths;
        results.differential_stress = results.sigma_t - results.sigma_r; % differential stress
        ds_max_depth = zeros(1,isave-1);
        for i=1:isave-1
            [~,ind] = max(results.differential_stress(:,i));
            ds_max_depth(i) = save_depths(ind);
        end
        figure, plot(results.time(mask)/seconds_in_year/1e6,ds_max_depth);
        ylabel('Depth of max. differential stress')

        %% plot a coulomb failure criterion
        tau_m = 0.5*abs(results.sigma_t-results.sigma_r);
        plith = zeros(size(results.sigma_t));
        for i=2:length(save_depths)
            if save_depths(i) <= h_crust
                rho1 = rhoc;
            else
                rho1 = rho;
            end
            plith(i,:) = plith(i-1,:) - rho1*g*(save_depths(i)-save_depths(i-1));
        end
        % plot coulomb failure criterion
        sigma_m = 0.5*(results.sigma_t+plith + results.sigma_r+plith);
        phi = atand(0.6);
        cohesion = 0.0;
        strength = cohesion*cosd(phi) - sigma_m*sind(phi);
        % Plot strength envelope from Mueller and Phillips 1995
        delta_sigma = -0.786*plith;


        xscale = 'linear';
        ax=[];
        figure();
        t=tiledlayout(7,1,'TileSpacing','compact','Padding','none');
        % t.Units = 'centimeters';
        % t.OuterPosition = [1 1 11 14];
        nexttile
        contourf(results.time(mask)/seconds_in_year/1e6,save_depths/1000,results.differential_stress(:,mask)/1e6,64,'Color','none'); %shading flat;
        hold on
        contour(results.time(mask)/seconds_in_year/1e6,save_depths/1000,results.differential_stress(:,mask)/1e6,[0 0],'k--'); %
        % contour(results.time(mask)/seconds_in_year/1e6,save_depths/1000,tau_m(:,mask) - strength(:,mask),[0 0],'Color','r','LineStyle','--'); %
        % contour(results.time(mask)/seconds_in_year/1e6,save_depths/1000,abs(results.differential_stress(:,mask))-delta_sigma(:,mask),[0 0],'Color','g','LineStyle','--'); %
        % contour(results.time(mask)/seconds_in_year/1e6,save_depths/1000,results.T(:,mask),[1000 1000],'Color','k','LineStyle','-'); %

        title(sprintf('Mars: eta0 %.2e, heating %.2f',results.state.mub,results.state.crust_heat_fraction))

        plot(results.time(mask)/seconds_in_year/1e6,((Ro-results.Ri(mask))+results.z(mask))/1000,'Color','k','LineWidth',1);
        %         set(gca,'YLim',[0 ceil(1+max(((Ro-results.Ri(mask))+results.z(mask))/1000))]);
        set(gca,'YDir','reverse');
        hcb = colorbar();
        set(gca,'Colormap',crameri('-roma'))
        stmax = max(max(abs(results.differential_stress(:,mask)/1e6)));
        caxis([-1 1]*stmax)
        hcb.Label.String = '\sigma_t-\sigma_r (MPa)';
        text(0.025,0.85,char('A'),'FontSize',12,'Units','normalized');
        % xlabel('Time (years)');
        % title(label);
        ylabel('Depth (km)');
        set(gca,'XScale',xscale);
        set(gca,'YLim',[0 50]);

        nexttile
        contourf(results.time(mask)/seconds_in_year/1e6,save_depths/1000,results.differential_stress(:,mask)/1e6,64,'Color','none'); %shading flat;

        contourf(results.time(mask)/seconds_in_year/1e6,save_depths/1000,results.sigma_t(:,mask)/1e6,64,'Color','none'); %shading flat;
        hold on
        plot(results.time(mask)/seconds_in_year/1e6,((Ro-results.Ri(mask))+results.z(mask))/1000,'Color','k','LineWidth',1);
        hold on
        contour(results.time(mask)/seconds_in_year/1e6,save_depths/1000,results.T(:,mask),[1000 1000],'Color','k','LineStyle','-'); %


        %         set(gca,'YLim',[0 ceil(1+max(((Ro-results.Ri(mask))+results.z(mask))/1000))]);
        set(gca,'YDir','reverse');
        hcb = colorbar();
        set(gca,'Colormap',crameri('-roma'))
        stmax = max(max(abs(results.sigma_t(:,mask)/1e6)));
        caxis([-1 1]*stmax)
        hcb.Label.String = '\sigma_t-\sigma_r (MPa)';
        text(0.025,0.85,char('A'+1),'FontSize',12,'Units','normalized');
        % xlabel('Time (years)');
        % title(label);
        ylabel('Depth (km)');
        set(gca,'XScale',xscale);
        hold on;
        % for i=1:ifail-1
        %     plot(results.failure_time(i)*1e6*[1 1],[results.failure_top(i) results.failure_bottom(i)]/1e3,'r');
        % end
        % TEMPERATURE
        nexttile
        contourf(results.time(mask)/seconds_in_year/1e6,save_depths/1000,results.T(:,mask),64,'Color','none'); %shading flat;
        hold on
        plot(results.time(mask)/seconds_in_year/1e6,((Ro-results.Ri(mask))+results.z(mask))/1000,'Color','k','LineWidth',1);
        %         set(gca,'YLim',[0 ceil(1+max(((Ro-results.Ri(mask))+results.z(mask))/1000))]);
        set(gca,'YDir','reverse');
        % ax1 = gca();
        % ax1.FontSize=8;
        set(gca,'Colormap',crameri('-lajolla'))
        hcb = colorbar();
        hcb.Label.String = 'Temperature (K)';
        text(0.025,0.85,char('C'),'FontSize',12,'Units','normalized');
        % xlabel('Time (years)');
        ylabel('Depth (km)');
        set(gca,'XScale',xscale);
        hold on;

        nexttile
        plot(results.time(mask)/seconds_in_year/1e6,results.ur(1,mask),'k')
        ylabel('u_r (m)')
        text(0.025,0.85,char('D'),'FontSize',12,'Units','normalized');

        nexttile
        plot(results.time(mask)/seconds_in_year/1e6,results.Pex(mask)/1e6,'k');
        ylabel('P_{ex} (MPa)');
        set(gca,'XScale',xscale);
        % ax2 = gca();
        % ax2.Position(3) = ax1.Position(3);
        % ax2.XLim = ax1.XLim;
        % ax2.FontSize=8;
        hold on
        % plot(results.failure_time(1:ifail-1)*1e6,results.failure_P(1:ifail-1)/1e6,'r.');
        % end_color = [0 0.9 0];
        % plot(results.failure_time(1:ifail-1)*1e6,(results.failure_P(1:ifail-1)+results.failure_dP(1:ifail-1))/1e6,'LineStyle','none','Color',end_color,'Marker','o','MarkerFaceColor',end_color,'MarkerSize',2);
        text(0.025,0.85,char('E'),'FontSize',12,'Units','normalized');
        % plot(results.time(mask)/seconds_in_year,results.Pex_crit(mask)/1e6,'k-');

        % xlabel('Time (years)');
        nexttile
        % for i=1:ifail-1
        %     if isnan(results.failure_erupted_volume(i))
        %         % plot nothing
        %     else
        %         if results.failure_P(i) - results.failure_Pex_crit(i) > 0
        %             plot(results.failure_time(i)*1e6*[1 1],[0 1],'b');
        %         else
        %             plot(results.failure_time(i)*1e6*[1 1],[0 1],'b--');
        %         end
        %     end
        % end
        plot(results.time(mask)/seconds_in_year/1e6,results.Tm(mask),'k-');
        text(0.025,0.85,char('F'),'FontSize',12,'Units','normalized');

        ylabel('T_m (K)');
        nexttile
        mask1 = mask & results.time>no_stress_time;
        plot(results.time(mask1)/seconds_in_year/1e6,results.stresss_crossover_depth(mask1)/1e3);
        hold on
        plot(results.time(mask1)/seconds_in_year/1e6,results.maximum_stress_depth(mask1)/1e3);
        legend('zero stress','maximum tension')
        ylabel('Depth (km)')
        text(0.025,0.85,char('G'),'FontSize',12,'Units','normalized');




        xlabel('Time (years)');
        set(gca,'XScale',xscale);
        % ax3=gca();
        % ax3.XLim = ax1.XLim;
        % ax3.Position(3) = ax1.Position(3);
        % ax3.Box = 'on';
        % ax3.FontSize=8;

        % set(gca,'XLim',[0 4500]);
        set(gca,'XLim',[0 4500]);

        fig = gcf();
        fig.Position(3:4) = [385   650];
        axmask = arrayfun(@(x) isa(x,'matlab.graphics.axis.Axes'),t.Children);
        linkaxes(t.Children(axmask),'x');
        

        set(t.Children(axmask),'XTickLabel',[])
        set(gca,'XTickLabel',get(gca,'XTick'))

        fig.Color = 'w';
        filename = sprintf('mars-thermal-evolution-zerotime-%f.pdf',no_stress_time/seconds_in_year/1e9);

    end
end