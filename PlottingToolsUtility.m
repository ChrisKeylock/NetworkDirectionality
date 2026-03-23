function PlottingToolsUtility(OutputFromMainCode)

% This plot is similar to panels (b), (c) and (d) of Fig. 1 in:
% Keylock and Carbone (2026), Phys. Rev. E 
% The input is the variable "output" from NetworkLaplacianDirectionality.m
% 


Ztau=OutputFromMainCode.Ztau;
d_eps=OutputFromMainCode.d_eps;

figure
subplot(2,2,[1 2])
eigDGL=eig(OutputFromMainCode.dgl);
plot(real(eigDGL),imag(eigDGL),'xc')
hold on
plot(real(Ztau{1}(:,1)),imag(Ztau{1}(:,1)),'-b','linewidth',0.7)
plot(real(Ztau{2}(:,1)),imag(Ztau{2}(:,1)),'--r','linewidth',0.7)
ylim([-0.01 max(imag(Ztau{1}(:,1)))*1.1])
xlabel('Real')
ylabel('imag')
legend('Eigenvalues', 'Z_{\tau}(L)', 'Z_{\tau}(L_{n})', 'Location', 'northeast')
str=strcat('\tau_{\epsilon} = ',num2str(OutputFromMainCode.epstau));
annotation('textbox',[0.13 0.83 0.2 0.08],'String',str,'FitBoxToText','on','LineStyle','none')

subplot(2,2,3)
plot(d_eps(:,2),d_eps(:,1),'-b','linewidth',0.7)
set(gca,'XDir','reverse')
xlabel('\theta_{0}')
ylabel('d_{\epsilon}')

subplot(2,2,4)
histobj=histogram(d_eps(:,1),round(sqrt(length(d_eps))),'normalization','probability')
ylabel('p(d_{\epsilon})')
xlabel('d_{\epsilon}')
ylim([0 1.05*max(histobj.Values)])
med_deps=round(median(d_eps(:,1)),4);
str2=strcat('< d_{\epsilon} >_{50} = ',num2str(med_deps));
annotation('textbox',[0.75 0.37 0.2 0.08],'String',str2,'FitBoxToText','on','LineStyle','none')






