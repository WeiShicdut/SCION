
figure


%%% CO2 (ppm)
subplot(4,3,1)
hold on
box on
plot(state.time_myr,state.CO2ppm,'k')
ylabel('Atm. CO_{2} (ppm)')

%%%% Ocean DIC
subplot(4,3,2)
hold on
box on
plot(state.time_myr,state.DIC_conc_p,'m')
plot(state.time_myr,state.DIC_conc_di,'k')
plot(state.time_myr,state.DIC_conc_s,'r')
plot(state.time_myr,state.DIC_conc_h,'c')
plot(state.time_myr,state.DIC_conc_d,'b')
legend('sm','dm','s','h','d','Location','southeast')
ylabel('DIC (mM)')

%%%% Ocean ALK
subplot(4,3,3)
hold on
box on
plot(state.time_myr,state.ALK_conc_p,'m')
plot(state.time_myr,state.ALK_conc_di,'k')
plot(state.time_myr,state.ALK_conc_s,'r')
plot(state.time_myr,state.ALK_conc_h,'c')
plot(state.time_myr,state.ALK_conc_d,'b')
ylabel('ALK (mM)')

%%%% Ocean pH
subplot(4,3,4)
hold on
box on
plot(state.time_myr,state.pH_p,'m')
plot(state.time_myr,state.pH_di,'k')
plot(state.time_myr,state.pH_s,'r')
plot(state.time_myr,state.pH_h,'c')
plot(state.time_myr,state.pH_d,'b')
ylabel('pH')


%%%% Ocean HCO3
subplot(4,3,5)
hold on
box on
plot(state.time_myr,state.HCO3_p,'m')
plot(state.time_myr,state.HCO3_di,'k')
plot(state.time_myr,state.HCO3_s,'r')
plot(state.time_myr,state.HCO3_h,'c')
plot(state.time_myr,state.HCO3_d,'b')
ylabel('HCO_{3}')

%%%% Ocean CO3
subplot(4,3,6)
hold on
box on
plot(state.time_myr,state.CO3_p,'m')
plot(state.time_myr,state.CO3_di,'k')
plot(state.time_myr,state.CO3_s,'r')
plot(state.time_myr,state.CO3_h,'c')
plot(state.time_myr,state.CO3_d,'b')
ylabel('CO_{3}')



%%%% Temperature
subplot(4,3,7)
hold on
box on
plot(state.time_myr,state.T_s,'r')
plot(state.time_myr,state.T_h,'c')
plot(state.time_myr,state.T_d,'b')
% plot(state.time_myr,state.T_cont,'c')
% plot(state.time_myr,state.GAST,'b')
legend('s','h','d','Location','southeast')
ylabel('Temperature (^oC)')

%%%% Fluxes
subplot(4,3,8)
hold on
box on
plot(state.time_myr,state.ccdeg,'m')
plot(state.time_myr,state.basw,'k')
plot(state.time_myr,state.granw,'r')
plot(state.time_myr,state.silw,'c')
plot(state.time_myr,state.carbw,'b')
legend('ccdeg','basw','granw','silw','carbw','Location','southeast')
ylabel('Fluxes (mol/yr)')

%%%% organic C burial
subplot(4,3,9)
hold on
box on
plot(state.time_myr,state.mocb_p,'m')
plot(state.time_myr,state.mocb_di,'k')
plot(state.time_myr,state.mocb_d,'b')
legend('p','di','d','Location','southeast')
ylabel('Organic C burial (mol/yr)')

%%%% Carbonate burial
subplot(4,3,10)
hold on
box on
plot(state.time_myr,state.mccb_p,'m')
plot(state.time_myr,state.mccb_di,'k')
plot(state.time_myr,state.mccb_d,'b')
legend('p','di','d','Location','southeast')
ylabel('Carbonate C burial (mol/yr)')
% 








% 
% % %%%% POC
% % subplot(4,6,9)
% % hold on
% % box on
% % plot(state.time_myr,state.POC_p,'m')
% % plot(state.time_myr,state.POC_di,'k')
% % plot(state.time_myr,state.POC_s,'r')
% % plot(state.time_myr,state.POC_h,'c')
% % plot(state.time_myr,state.POC_d,'b')
% % %legend('sm','dm','s','h','d','Location','southeast')
% % ylabel('POC (mol)')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%% DP
% subplot(4,6,10)
% hold on
% box on
% plot(state.time_myr,log10(state.DP_conc_p),'m')
% plot(state.time_myr,log10(state.DP_conc_di),'k')
% plot(state.time_myr,log10(state.DP_conc_s),'r')
% plot(state.time_myr,log10(state.DP_conc_h),'c')
% plot(state.time_myr,log10(state.DP_conc_d),'b')
% %legend('sm','dm','s','h','d','Location','southeast')
% ylabel('log(DP) (mM)')
% 
% %%%% O2
% subplot(4,6,11)
% hold on
% box on
% plot(state.time_myr,log10(state.O2_conc_p),'m')
% plot(state.time_myr,log10(state.O2_conc_di),'k')
% plot(state.time_myr,log10(state.O2_conc_s),'r')
% plot(state.time_myr,log10(state.O2_conc_h),'c')
% plot(state.time_myr,log10(state.O2_conc_d),'b')
% %legend('sm','dm','s','h','d','Location','southeast')
% ylabel('log(O_2) conc (mM)')
% 
% 
% %%%% FeIII
% subplot(4,6,12)
% hold on
% box on
% plot(state.time_myr,log10(state.FeIII_conc_pf),'m')
% plot(state.time_myr,log10(state.FeIII_conc_dif),'k')
% plot(state.time_myr,log10(state.FeIII_conc_sf),'r')
% plot(state.time_myr,log10(state.FeIII_conc_hf),'c')
% plot(state.time_myr,log10(state.FeIII_conc_df),'b')
% %legend('sm','dm','s','h','d','Location','southeast')
% xlabel(str)
% ylabel('log(FeIII) conc (mM)')
% 
% %%%% SO4
% subplot(4,6,13)
% hold on
% box on
% plot(state.time_myr,state.SO4_conc_pf,'m')
% plot(state.time_myr,state.SO4_conc_dif,'k')
% plot(state.time_myr,state.SO4_conc_sf,'r')
% plot(state.time_myr,state.SO4_conc_hf,'c')
% plot(state.time_myr,state.SO4_conc_df,'b')
% %legend('sm','dm','s','h','d','Location','southeast')
% xlabel(str)
% ylabel('SO_4 conc (mM)')
% 
% %%%% FeII
% subplot(4,6,14)
% hold on
% box on
% plot(state.time_myr,log10(state.FeII_conc_pf),'m')
% plot(state.time_myr,log10(state.FeII_conc_dif),'k')
% plot(state.time_myr,log10(state.FeII_conc_sf),'r')
% plot(state.time_myr,log10(state.FeII_conc_hf),'c')
% plot(state.time_myr,log10(state.FeII_conc_df),'b')
% %legend('sm','dm','s','h','d','Location','southeast')
% xlabel(str)
% ylabel('log(FeII) conc (mM)')
% 
% %%%% H2S
% subplot(4,6,15)
% hold on
% box on
% plot(state.time_myr,log10(state.H2S_conc_pf)+3,'m')
% plot(state.time_myr,log10(state.H2S_conc_dif)+3,'k')
% plot(state.time_myr,log10(state.H2S_conc_sf)+3,'r')
% plot(state.time_myr,log10(state.H2S_conc_hf)+3,'c')
% plot(state.time_myr,log10(state.H2S_conc_df)+3,'b')
% %legend('sm','dm','s','h','d','Location','southeast')
% xlabel(str)
% ylabel('log(H2S) conc (\muM)')
% 
% %%%% FeIIIw
% subplot(4,6,16)
% hold on
% box on
% plot(state.time_myr,state.FeIIIwf,'m')
% plot(state.time_myr,state.phoswf,'b')
% legend('FeIII weathering','phosw','Location','southeast')
% xlabel(str)
% ylabel('FWeathering fluxes (mol/yr)')
% 
% %%%% FeIIIscavenging
% subplot(4,6,17)
% hold on
% box on
% plot(state.time_myr,state.FeIIIscavenging_pf,'m')
% plot(state.time_myr,state.FeIIIscavenging_dif,'k')
% plot(state.time_myr,state.FeIIIscavenging_sf,'r')
% plot(state.time_myr,state.FeIIIscavenging_hf,'c')
% plot(state.time_myr,state.FeIIIscavenging_df,'b')
% %legend('sm','dm','s','h','d','Location','southeast')
% xlabel(str)
% ylabel('FeIII scavenging (mol/yr)')
% 
% %%%% S cycle
% subplot(4,6,18)
% hold on
% box on
% plot(state.time_myr,state.pyritewf,'m')
% plot(state.time_myr,state.sulfatewf,'k')
% %plot(state.time_myr,state.pyriteb_pf,'r')
% %plot(state.time_myr,state.pyriteb_dif,'c')
% plot(state.time_myr,state.sulfatebf,'b')
% legend('pyritew','sulfatew','sulfateb','Location','southeast')
% xlabel(str)
% ylabel('S fluxes (mol/yr)')
% 
% %%%% pyrite formation
% subplot(4,6,19)
% hold on
% box on
% plot(state.time_myr,state.pyF_pf,'m')
% plot(state.time_myr,state.pyF_dif,'k')
% plot(state.time_myr,state.pyF_sf,'r')
% plot(state.time_myr,state.pyF_hf,'c')
% plot(state.time_myr,state.pyF_df,'b')
% %legend('sm','dm','s','h','d','Location','southeast')
% xlabel(str)
% ylabel('Pyrite formation (mol/yr)')
% 
% %%%% iron oxidation
% subplot(4,6,20)
% hold on
% box on
% plot(state.time_myr,state.ironO_pf,'m')
% plot(state.time_myr,state.ironO_dif,'k')
% plot(state.time_myr,state.ironO_sf,'r')
% plot(state.time_myr,state.ironO_hf,'c')
% plot(state.time_myr,state.ironO_df,'b')
% %legend('sm','dm','s','h','d','Location','southeast')
% xlabel(str)
% ylabel('Iron oxidation (mol/yr)')
% 
% %%%% iron reduction with sulfur oxidation
% subplot(4,6,21)
% hold on
% box on
% plot(state.time_myr,state.SironR_pf,'m')
% plot(state.time_myr,state.SironR_dif,'k')
% plot(state.time_myr,state.SironR_sf,'r')
% plot(state.time_myr,state.SironR_hf,'c')
% plot(state.time_myr,state.SironR_df,'b')
% %legend('sm','dm','s','h','d','Location','southeast')
% xlabel(str)
% ylabel('Iron red S oxi (mol/yr)')
% 
% %%%% SideP
% subplot(4,6,22)
% hold on
% box on
% plot(state.time_myr,state.mocb_p_FeIIIf,'m')
% plot(state.time_myr,state.ocbother_pf,'c')
% legend('mocbFeIII_p','ocbscavenging_p','ocbother_p','Location','southeast')
% xlabel(str)
% ylabel('OC burial (mol/yr)')
% 
% %%%% pripr
% subplot(4,6,23)
% hold on
% box on
% plot(state.time_myr,state.mocb_di_FeIIIf,'m')
% plot(state.time_myr,state.ocbother_dif,'c')
% legend('mocbFeIII_d_i','ocbscavenging_d_i','ocbother_d_i','Location','southeast')
% xlabel(str)
% ylabel('OC burial (mol/yr)')
% 
% %%%% remin
% subplot(4,6,24)
% hold on
% box on
% plot(state.time_myr,state.mocb_d_FeIIIf,'m')
% plot(state.time_myr,state.ocbother_df,'c')
% legend('mocbFeIII_d','ocbscavenging_d','ocbother_d','Location','southeast')
% xlabel(str)
% ylabel('OC burial (mol/yr)')

