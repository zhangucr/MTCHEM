% MTCHEM.m
% alpha-pinene ozonolysis mechanisms
% Authors: Haofei Zhang and Chuanyang Shen (UC Riverside)
% Last update: July 2025

% General tunable parameters
fautox = 0.1; % tuning parameter for autoxidaB tion rate constant scaling. suggest: 0.01 - 1.
fROOH_sec = 1; % fraction of ROOH formation from secondary RO2 + HO2.
fROOH_tert = 0.2; % fraction of ROOH formation from tertiary RO2 + HO2. AcylRO2 is not included.
fsec = 0.5; % fraction of a RO2 being secondary. Now only used for APC10H15O8
ftert = 0.5; % fraction of a RO2 being tertiary. Now only used for APC10H15O8
facyl = 1-fsec-ftert; % fraction of a RO2 being acyl. Now only used for APC10H15O8
fdecomp_C109O = 0.9; % fraction of C109O decomposition
fdecomp_C89CO2 = 0.6; % fraction of C89CO2 decomposition
KRO2RO2_fast = 2E-12;% can be as high as 1e-10, based on Berdnt et al., 2018.
fAPINOOB_RB = 0.1; %branching ratio of APINOOdecomposition into ring-broken C10H15O4-RO2. Yield from a-pinene +O3 ranges from 0.03-0.25 based on Iyer et al., 2021 Nature Comm, translated to 0.1 - 0.66 for this factor

kautox_gen = 0.1.*fautox; % generic RO2 autoxidation rate constant, s-1.
kselfterm = 0.1; % RO2 self-termination rate constant, s-1.
bRO_ps = 0.6; % branching ratio of RO formation for primary and secondary RO2 + RO2 reactions. Default value is from MCM but can go up to 0.8
bROH_ps = 0.2; % branching ratio of ROH (and R=O) formation for primary and secondary RO2 + RO2 reactions. Default value is from MCM.
bRO_ta = 0.7; % branching ratio of RO formation for tertiary and acyl RO2 + RO2 reactions. Default value is from MCM but can go up to 0.8
bROH_ta = 0.3; % branching ratio of ROH formation for tertiary and acyl RO2 + RO2 reactions.Default value is from MCM.

% Autoxidation rate constants reported in literature:
kautox_PINAL4RO2 = 1.08E-66.*T^25.23.*exp(1616./T).*fautox; % 1,4-aldehydic H shift. rate constant based on SAR. k298 ~ 0.065
kautox_C109O2 = 3.89E-4.*T^4.39.*exp(-6565./T).*fautox; % 1,5-H shift. rate constant based on SAR. k298~0.008
kautox_C107O2 = 1.46E-30.*T^12.9.*exp(-2663./T).*fautox; % 1,7-aldehydic H shift. rate constant based on SAR. k298 ~ 0.016. Assuming 50% of C107O2 is in cis-configuration that can undergo this isomerization
kautox_APC10H15O4RB = 20.*fautox; % rate constant based on Iyer et al., 2020 Nature Comm. Using a lower value here consistent with Zhao et al., 2021 ES&T.
kautox_PINAL4H1RO2 = 4.03E-2.*T^3.81.*exp(-6030./T).*fautox; % based on SAR (Vereecken and Noziere, 2020), 1,6-H (a-OH) shift. Assuming 50% of this RO2 is in cis configuration. k298~0.175
kautox_C108O2 = 2.03E-42.*T^16.79.*exp(-268./T).*fautox; % based on SAR (Vereecken and Noziere, 2020), 1,6-H (aldehydic) shift.k298 ~0.28
kautox_C920CO3 = 4.48E-15.*T^8.73.*exp(-6280./T).*fautox; % based on SAR (Vereecken and Noziere, 2020), acylRO2 correction, 1,5-H shift. k298~0.013
kautox_PINALPA4RO2 = 6.40E-35.*T^14.71.*exp(-1831./T).*fautox; % based on SAR (Vereecken and Noziere, 2020), 1,6-H shift with exo-beta-oxo correction. This RO2 is formed from rapid scrambling of PINAL4P3RO2. k298~0.34. Assuming 50% cis.
kautox_PINALPA1RO2 = 2.24E-16.*T^8.73.*exp(-6280./T).*fautox; % based on SAR (Vereecken and Noziere, 2020), 1,5-H shift. Slow. This RO2 is formed from rapid scrambling of PINAL1P3RO2. k298~0.001
kautox_PINAL10P1RO2 = 1.46E-30.*T^12.9.*exp(-2663./T).*fautox; % 1,7-aldehydic H shift. rate constant based on SAR. k298 ~ 0.016. Assuming 50% of C107O2 is in cis-configuration that can undergo this isomerization
kautox_APC10H15O6RO2 = 1.28E-34.*T^14.71.*exp(-1831./T).*fautox; % based on SAR (Vereecken and Noziere, 2020),1,6-H shift with exo-beta-oxo correction. k298 ~ 0.68
kautox_APC10H15O7RO2 = 3.50E-21.*T^10.79.*exp(-5440./T).*fautox; % based on SAR (Vereecken and Noziere, 2020),dominated by tertRO2 1,5-H shift with exo-beta-oxo. k298 ~ 0.021
kautox_APC10H15O8RO2 = 7.56E-67.*T^25.23.*exp(1616./T).*fautox; % based on SAR (Vereecken and Noziere, 2020), 1,4-aldehydic H shift. ~ 70% isomers are this type of RO2. The others don't autoxidize. k298~0.05
kautox_APC10H15O10RO2 = 0.0.*fautox; % structure dependent. highly uncertainty
kautox_C10H15O3O2 = 3.0.*fautox; % from PRAM, 298K
kautox_C10H15O4O2 = 3.0.*fautox; % from PRAM, 298K
kautox_C10H15O5O2 = 1.5.*fautox; % from PRAM, 298K
kautox_C10H15O6O2 = 0.5.*fautox; % from PRAM, 298K

% The following five rate constants are for RO2+RO2=ROOR reactions. Five
% rate constants are used for different level of oxidized RO2:
kROORa = 1e-14; % use kROORa if product ROOR has nO 2-6
kROORb = 5e-14; % use kROORb if product ROOR has nO 7-9
kROORc = 1e-13; % use kROORb if product ROOR has nO 10-12
kROORd = 5e-13; % use kROORb if product ROOR has nO 13-15
kROORe = 1e-12; % use kROORb if product ROOR has nO 16-18
fROOR = 1.0; % tuning factor for kROOR
fROOR1 = 1.0; % tuning factor for BU2TOL dimer formation
fROOR2 = 1.0; % tuning factor for CHEX dimer formation

SpeciesToAdd = {...
% Criegee intermediate related products:
'APINAOO'; % SCI from APINOOA;
'AAHP';% lumped AAHP products from SCI + acids/alcohols.
'SOZ';% SOZ from SCI + acetone.
% C10-RO2 species:
'PINAL4RO2'; % C10H15O4
'APC10H15O4RB'; % ring-opened C10H15O4
'PINAL4H1RO2'; % C10H15O5
'APC10H15O5RO2'; % C10H15O5
'PINALPA4RO2'; % C10H15O6. Formed from PINAL4P3RO2 rapid scrambling.
'APC10H15O6RO2'; % C10H15O6, mainly from CAPC10H15O4RB autoxidation
'PINAL10P1RO2'; % C10H15O6
'PINALPA1RO2'; % C10H15O6. Formed from PINAL1P3RO2 rapid scrambling.
'APC10H15O7RO2';
'APC10H15O8RO2';
'APC10H15O9RO2';
'APC10H15O10RO2';
% C10-RO species:
'PINAL4RO';
'PINAL4H1RO';
'PINALPA4RO';
'PINALPA1RO';
'PINAL10P1RO';
'APC10H15O2RO';
'APC10H15O3RO';
'APC10H15O4RO';
'APC10H15O5RO';
'APC10H15O6RO';
'APC10H15O7RO';
'APC10H15O8RO';
'APC10H15O9RO';
% C10-CHO Closed-shell species:
'PINAL4P'; % C10H16O4; 
'PINAL4H'; % C10H16O3; 
'PINAL4C'; % C10H14O3; 
'PINAL4H1P'; % C10H16O5; 
'PINAL1H4H'; % C10H16O4; 
'APC10H14O'; % 
'APC10H14O2'; % 
'APC10H14O4'; % 
'APC10H14O5'; % 
'APC10H14O6'; % 
'APC10H14O7'; % 
'APC10H14O8'; % 
'APC10H14O9'; % 
'APC10H16O'; % 
'APC10H16O2'; % 
'APC10H16O3'; % 
'APC10H16O4'; % 
'APC10H16O5'; % 
'APC10H16O6'; % 
'APC10H16O7'; % 
'APC10H16O8'; % 
'APC10H16O9'; % 
'APC10H16O10'; % 
% C10-N-containing species:
'APC10H15NO5'; % 
'APC10H15NO6'; % 
'APC10H15NO7'; % 
'APC10H15NO8'; % 
'APC10H15NO9'; % 
'APC10H15NO10'; % 
'APC10H15NO11'; % 

% dimer species:
'C20H30O2';'C20H30O3';'C20H30O4';'C20H30O5';'C20H30O6';'C20H30O7';'C20H30O8';'C20H30O9';'C20H30O10';'C20H30O11';'C20H30O12';'C20H30O13';'C20H30O14';'C20H30O15';'C20H30O16';'C20H30O17';'C20H30O18';...
'C20H32O3';'C20H32O4';'C20H32O5';'C20H32O6';'C20H32O7';'C20H32O8';'C20H32O9';'C20H32O10';'C20H32O11';'C20H32O12';'C20H32O13';'C20H32O14';'C20H32O15';'C20H32O16';'C20H32O17';'C20H32O18';...
'C20H34O4';'C20H34O5';'C20H34O6';'C20H34O7';'C20H34O8';'C20H34O9';'C20H34O10';'C20H34O11';'C20H34O12';'C20H34O13';'C20H34O14';'C20H34O15';'C20H34O16';'C20H34O17';'C20H34O18';...
'C20H29O3NO3';'C20H29O4NO3';'C20H29O5NO3';'C20H29O6NO3';'C20H29O7NO3';'C20H29O8NO3';'C20H29O9NO3';'C20H29O10NO3';'C20H29O11NO3';'C20H29O12NO3';...
'C20H31O2NO3';'C20H31O3NO3';'C20H31O4NO3';'C20H31O5NO3';'C20H31O6NO3';'C20H31O7NO3';'C20H31O8NO3';'C20H31O9NO3';'C20H31O10NO3';'C20H31O11NO3';'C20H31O12NO3';'C20H31O13NO3';'C20H31O14NO3';'C20H31O15NO3';...
'C20H33O4NO3';'C20H33O5NO3';'C20H33O6NO3';'C20H33O7NO3';'C20H33O8NO3';'C20H33O9NO3';'C20H33O10NO3';'C20H33O11NO3';'C20H33O12NO3';...
'C20H35O4NO3';'C20H35O5NO3';'C20H35O6NO3';'C20H35O7NO3';'C20H35O8NO3';'C20H35O9NO3';'C20H35O10NO3';'C20H35O11NO3';...
'C19H26O4';'C19H26O5';'C19H26O6';'C19H26O7';'C19H26O8';'C19H26O9';'C19H26O10';'C19H26O11';'C19H26O12';'C19H26O13';'C19H26O14';...
'C19H28O3';'C19H28O4';'C19H28O5';'C19H28O6';'C19H28O7';'C19H28O8';'C19H28O9';'C19H28O10';'C19H28O11';'C19H28O12';'C19H28O13';'C19H28O14';'C19H28O15';'C19H28O16';'C19H28O17';'C19H28O18';...
'C19H30O3';'C19H30O4';'C19H30O5';'C19H30O6';'C19H30O7';'C19H30O8';'C19H30O9';'C19H30O10';'C19H30O11';'C19H30O12';'C19H30O13';'C19H30O14';...
'C19H32O4';'C19H32O5';'C19H32O6';'C19H32O7';'C19H32O8';'C19H32O9';'C19H32O10';'C19H32O11';'C19H32O12';'C19H32O13';'C19H32O14';...
'C19H29O3NO3';'C19H29O4NO3';'C19H29O5NO3';'C19H29O6NO3';'C19H29O7NO3';'C19H29O8NO3';'C19H29O9NO3';'C19H29O10NO3';'C19H29O11NO3';'C19H29O12NO3';'C19H29O13NO3';'C19H29O14NO3';'C19H29O15NO3';...
'C19H31O5NO3';'C19H31O6NO3';'C19H31O7NO3';'C19H31O8NO3';'C19H31O9NO3';'C19H31O10NO3';'C19H31O11NO3';...
'C19H33O5NO3';'C19H33O6NO3';'C19H33O7NO3';'C19H33O8NO3';'C19H33O9NO3';'C19H33O10NO3';'C19H33O11NO3';...
'C18H26O2';'C18H26O3';'C18H26O4';'C18H26O5';'C18H26O6';'C18H26O7';'C18H26O8';'C18H26O9';'C18H26O10';'C18H26O11';'C18H26O12';'C18H26O13';'C18H26O14';'C18H26O15';'C18H26O16';'C18H26O17';'C18H26O18';...
'C18H28O2';'C18H28O3';'C18H28O4';'C18H28O5';'C18H28O6';'C18H28O7';'C18H28O8';'C18H28O9';'C18H28O10';'C18H28O11';'C18H28O12';'C18H28O13';'C18H28O14';...
'C18H30O4';'C18H30O5';'C18H30O6';'C18H30O7';'C18H30O8';'C18H30O9';'C18H30O10';'C18H30O11';'C18H30O12';'C18H30O13';'C18H30O14';...
'C18H27O4NO3';'C18H27O5NO3';'C18H27O6NO3';'C18H27O7NO3';'C18H27O8NO3';'C18H27O9NO3';'C18H27O10NO3';'C18H27O11NO3';'C18H27O12NO3';'C18H27O13NO3';'C18H27O14NO3';'C18H27O15NO3';'C18H27O16NO3';...
'C18H29O6NO3';'C18H29O7NO3';'C18H29O8NO3';'C18H29O9NO3';'C18H29O10NO3';'C18H29O11NO3';'C18H29O12NO3';...
'C18H31O6NO3';'C18H31O7NO3';'C18H31O8NO3';'C18H31O9NO3';'C18H31O10NO3';'C18H31O11NO3';'C18H31O12NO3';...
'C17H24O7';'C17H24O8';...
'C17H26O5';'C17H26O6';'C17H26O7';'C17H26O8';'C17H26O9';...
'C17H28O4';'C17H28O5';'C17H28O6';'C17H28O7';'C17H28O8';'C17H28O9';'C17H28O10';...
'C16H22O8';...
'C16H24O6';'C16H24O7';'C16H24O8';'C16H24O9';...
'C16H26O4';'C16H26O5';'C16H26O6';'C16H26O7';'C16H26O8';'C16H26O9';'C16H26O10';...
'BUT2OLdimer';'CHEXdimer';

% fragmentation products:
% RO2s:
'C7RO2';'C8RO2';'C9RO2';'fragRO2';'C7NRO2';'C8NRO2';'C9NRO2';'fragNRO2';
% closed-shell products:
'cycAPHA'; % the cyclic acylperoxyhemiacetal from C89CO3H, intermediate for cis-pinic acid
'C9OOH';'C9H';'C9C'; % Lumped C9-CHO compounds. Average formula = C9H14O6. Average MW = 218. Average vapor pressure (298K) = 9.61e-12 atm. Average vapor pressure (400K) = 1.42e-6 atm.
'C9NO3';'C9NOOH';'C9NH';'C9NC'; % Lumped C9-CHNO compounds. Average formula = C9H13.7NO7.6. Average MW = 257. Average vapor pressure (298K) = 1.61E-11 atm. Average vapor pressure (400K) = 1.53e-6 atm.
'C8OOH';'C8H';'C8C';% Lumped C8-CHO compounds. Average formula = C8H11.4O5.4. Average MW = 194. Average vapor pressure (298K) = 1.17e-10 atm. Average vapor pressure (400K) = 7.96e-6 atm.
'C8NO3';'C8NOOH';'C8NH';'C8NC';% Lumped C8-CHNO compounds. Average formula = C8H11.5NO7.2. Average MW = 236. Average vapor pressure (298K) = 1.21e-10 atm. Average vapor pressure (400K) = 7.04e-6 atm.
'C7OOH';'C7H';'C7C';% Lumped C7-CHO compounds. Average formula = C7H10.5O5.2. Average MW = 178. Average vapor pressure (298K) = 8.53e-10 atm. Average vapor pressure (400K) = 2.89e-5 atm.
'C7NO3';'C7NOOH';'C7NH';'C7NC';% Lumped C7-CHNO compounds. Average formula = C7H10.2NO7.4. Average MW = 227. Average vapor pressure (298K) = 1.91E-10 atm. Average vapor pressure (400K) = 8.78e-6 atm.
'fragOOH';'fragH';'fragC';% Lumped <C7 CHO compounds. Average formula = C6H9.3O5.6. Average MW = 170. Average vapor pressure (298K) = 3.03e-10 atm. Average vapor pressure (400K) = 1.57e-5 atm.
'fragNO3';'fragNOOH';'fragNH';'fragNC';% Lumped <C7 CHNO compounds. Average formula = C6H9.6NO7.4. Average MW = 214. Average vapor pressure (298K) = 5.07e-10 atm. Average vapor pressure (400K) = 1.69e-5 atm.
% process species as indicators of reaction types:
'autox';'term';'O2';'isom';'decomp';
};

RO2ToAdd = {...
'C107O2';'C109O2';...
'PINAL4RO2';'APC10H15O4RB';'PINALPA4RO2';'APC10H15O6RO2';...
'PINAL10P1RO2';'PINALPA1RO2';'PINAL4H1RO2';'C7RO2';'APC10H15O5RO2';'APC10H15O7RO2';'C8RO2';'APC10H15O8RO2';...
'C9RO2';'fragRO2';'APC10H15O9RO2';'APC10H15O10RO2';'APC10H15O9RO2';...
'C7NRO2';'C8NRO2';'C9NRO2';'fragNRO2';...
};

AddSpecies

%% Modified alpha-pinene oxidation.
% new branching ratios based on recent studies (Claflin et al., 2018 ACS Earth and Space Chem. and Berndt et al., 2016 Nature Comm.)

RxnToReplace = 'APINENE + O3 = APINOOA';
kToReplace = 8.05E-16.*exp(-640./T).*0.6;
ReplaceRxn

RxnToReplace = 'APINENE + O3 = APINOOB';
kToReplace = 8.05E-16.*exp(-640./T).*0.4;
ReplaceRxn

%% Formation of PINIC ACID precursor
% cycAPHA assumed to produce pinic-acid in the condensed phase
i=i+1;
Rnames{i} = 'C89CO3H = cycAPHA';% New reaction to represent intramolecular BV reaction to form pinic acid through a cyclic acylperoxyhemiacetal (Zhang and Zhang, 2021 Anal. Chem.)
k(:,i) = 1.0; % rate consant is unknown, but guessed to be 1 s^-1 by HZ
Gstr{i,1} = 'C89CO3H'; 
fC89CO3H(i)=fC89CO3H(i)-1; fcycAPHA(i)=fcycAPHA(i)+1;

i=i+1;
Rnames{i} = 'cycAPHA = C89CO3H'; % reverse reaction
k(:,i) = 0.1; % rate consant is unknown, but guessed to be 0.1 s^-1 by HZ. Modeled SOA not very sensitive to this value.
Gstr{i,1} = 'cycAPHA'; 
fcycAPHA(i)=fcycAPHA(i)-1; fC89CO3H(i)=fC89CO3H(i)+1; 
%% Modified RO2+HO2=ROOH yields
% this part of the mechanism can be applied to MCM
% The following changes from MCM consider that tertiary RO2 + HO2 do not
% form ROOH at 100% yield. The tertiary RO2 include C107O2, C108O2, C921O2,
% C97O2
RxnToReplace = 'C107O2 + HO2 = C107OOH';
kToReplace = KRO2HO2.*0.914.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'C107O2 + HO2 = C107O';
k(:,i) = KRO2HO2.*0.914.*(1-fROOH_tert);
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'HO2'; 
fC107O2(i)=fC107O2(i)-1; fHO2(i)=fHO2(i)-1; fC107O(i)=fC107O(i)+1;

RxnToReplace = 'C108O2 + HO2 = C108OOH';
kToReplace = KRO2HO2.*0.914.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'C108O2 + HO2 = C108O';
k(:,i) = KRO2HO2.*0.914.*(1-fROOH_tert);
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'HO2'; 
fC108O2(i)=fC108O2(i)-1; fHO2(i)=fHO2(i)-1; fC108O(i)=fC108O(i)+1;

RxnToReplace = 'C812O2 + HO2 = C812OOH';
kToReplace = KRO2HO2.*0.859.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'C812O2 + HO2 = C812O';
k(:,i) = KRO2HO2.*0.859.*(1-fROOH_tert);
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'HO2'; 
fC812O2(i)=fC812O2(i)-1; fHO2(i)=fHO2(i)-1; fC812O(i)=fC812O(i)+1; 

RxnToReplace = 'C921O2 + HO2 = C921OOH';
kToReplace = KRO2HO2.*0.890.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'C921O2 + HO2 = C921O';
k(:,i) = KRO2HO2.*0.890.*(1-fROOH_tert);
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'HO2'; 
fC921O2(i)=fC921O2(i)-1; fHO2(i)=fHO2(i)-1; fC921O(i)=fC921O(i)+1; 

RxnToReplace = 'C97O2 + HO2 = C97OOH';
kToReplace = KRO2HO2.*0.890.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'C97O2 + HO2 = C97O';
k(:,i) = KRO2HO2.*0.890.*(1-fROOH_tert);
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'HO2'; 
fC97O2(i)=fC97O2(i)-1; fHO2(i)=fHO2(i)-1; fC97O(i)=fC97O(i)+1;

RxnToReplace = 'PINALO2 + HO2 = PINALOOH';
kToReplace = KRO2HO2.*0.914.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'PINALO2 + HO2 = PINALO';
k(:,i) = KRO2HO2.*0.914.*(1-fROOH_tert);
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'HO2'; 
fPINALO2(i)=fPINALO2(i)-1; fHO2(i)=fHO2(i)-1; fPINALO(i)=fPINALO(i)+1; 

RxnToReplace = 'C106O2 + HO2 = C106OOH';
kToReplace = KRO2HO2.*0.914.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'C106O2 + HO2 = C106O';
k(:,i) = KRO2HO2.*0.914.*(1-fROOH_tert);
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'HO2'; 
fC106O2(i)=fC106O2(i)-1; fHO2(i)=fHO2(i)-1; fC106O(i)=fC106O(i)+1; 

RxnToReplace = 'APINAO2 + HO2 = APINAOOH';
kToReplace = KRO2HO2.*0.914.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'APINAO2 + HO2 = APINAO';
k(:,i) = KRO2HO2.*0.914.*(1-fROOH_tert);
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'HO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fHO2(i)=fHO2(i)-1; fAPINAO(i)=fAPINAO(i)+1;

RxnToReplace = 'APINBO2 + HO2 = APINBOOH';
kToReplace = KRO2HO2.*0.914.*fROOH_sec;
ReplaceRxn
i=i+1;
Rnames{i} = 'APINBO2 + HO2 = APINBO';
k(:,i) = KRO2HO2.*0.914.*(1-fROOH_sec);
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'HO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fHO2(i)=fHO2(i)-1; fAPINBO(i)=fAPINBO(i)+1;

RxnToReplace = 'APINCO2 + HO2 = APINCOOH';
kToReplace = KRO2HO2.*0.914.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'APINCO2 + HO2 = APINCO';
k(:,i) = KRO2HO2.*0.914.*(1-fROOH_tert);
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'HO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fHO2(i)=fHO2(i)-1; fAPINCO(i)=fAPINCO(i)+1; 

RxnToReplace = 'C98O2 + HO2 = C98OOH';
kToReplace = KRO2HO2.*0.890.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'C98O2 + HO2 = C98O';
k(:,i) = KRO2HO2.*0.890.*(1-fROOH_tert);
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'HO2'; 
fC98O2(i)=fC98O2(i)-1; fHO2(i)=fHO2(i)-1; fC98O(i)=fC98O(i)+1;

RxnToReplace = 'C922O2 + HO2 = C922OOH';
kToReplace = KRO2HO2.*0.890.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'C922O2 + HO2 = C922O';
k(:,i) = KRO2HO2.*0.890.*(1-fROOH_tert);
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'HO2'; 
fC922O2(i)=fC922O2(i)-1; fHO2(i)=fHO2(i)-1; fC922O(i)=fC922O(i)+1; 

RxnToReplace = 'C86O2 + HO2 = C86OOH';
kToReplace = KRO2HO2.*0.859.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'C86O2 + HO2 = C86O';
k(:,i) = KRO2HO2.*0.859.*(1-fROOH_tert);
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'HO2'; 
fC86O2(i)=fC86O2(i)-1; fHO2(i)=fHO2(i)-1; fC86O(i)=fC86O(i)+1;

RxnToReplace = 'C810O2 + HO2 = C810OOH';
kToReplace = KRO2HO2.*0.914.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'C810O2 + HO2 = C810O';
k(:,i) = KRO2HO2.*0.914.*(1-fROOH_tert);
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'HO2'; 
fC810O2(i)=fC810O2(i)-1; fHO2(i)=fHO2(i)-1; fC810O(i)=fC810O(i)+1;

RxnToReplace = 'C813O2 + HO2 = C813OOH';
kToReplace = KRO2HO2.*0.859.*fROOH_tert;
ReplaceRxn
i=i+1;
Rnames{i} = 'C813O2 + HO2 = C813O';
k(:,i) = KRO2HO2.*0.859.*fROOH_tert;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'HO2'; 
fC813O2(i)=fC813O2(i)-1; fHO2(i)=fHO2(i)-1; fC813O(i)=fC813O(i)+1;

%% Modified RO's fates for major MCM-ROs based on SAR
% C89CO2's fate is split at 20:80 to keep pinic acid amount ~ 5-10% of total
% SOA. MCM assumes 80:20 while SAR suggests 0:100.
RxnToReplace = 'C89CO2 = C811CO3';
kToReplace = KDEC.*(1-fdecomp_C89CO2);
ReplaceRxn
RxnToReplace = 'C89CO2 = C89O2';
kToReplace = KDEC.*fdecomp_C89CO2;
ReplaceRxn

% C109O fate is split at 100:0 for now (decomp. vs. isom.). MCM assumes
% 80:20, while SAR suggests 100:0.
RxnToReplace = 'C109O = C89CO3 + HCHO';
kToReplace = KDEC.*fdecomp_C109O;
ReplaceRxn
RxnToReplace = 'C109O = C920CO3'; % 1,5-H shift is more likely to happen, but still outrun by decomposition based on SAR
kToReplace = KDEC.*(1-fdecomp_C109O);
ReplaceRxn

RxnToReplace = 'C108O = C717O2 + CH3COCH3';
kToReplace = KDEC.*0.92;
ReplaceRxn
i=i+1;
Rnames{i} = 'C108O = APC10H15O6RO2';% 1,5-H shift (aldehydic H) to form a acylRO2. For now, use the generic formula. But could be explicitly represented later.
k(:,i) = KDEC.*0.08;
Gstr{i,1} = 'C108O'; 
fC108O(i)=fC108O(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)+1;

RxnToReplace = 'APINCO = CH3COCH3 + C720O2';
kToReplace = KDEC.*0.87;
ReplaceRxn
i=i+1;
Rnames{i} = 'APINCO = APC10H17O4RO2';%1,5-H shift (alpha-OH). But Claeys et al., 2009 EST suggested the other 1,5-H can be fast too and may lead to terpenylic acid
k(:,i) = KDEC.*0.13;
Gstr{i,1} = 'APINCO'; 
fAPINCO(i)=fAPINCO(i)-1; fAPC10H17O4RO2(i)=fAPC10H17O4RO2(i)+1;

RxnToReplace = 'C86O = C511O2 + CH3COCH3';
kToReplace = KDEC.*0.92;
ReplaceRxn
i=i+1;
Rnames{i} = 'C86O = C8RO2'; % 1,5-H shift (aldehydic H) to produce C8H13O5RO2. For now, use the generic formula for product.
k(:,i) = KDEC.*0.08;
Gstr{i,1} = 'C86O'; 
fC86O(i)=fC86O(i)-1; fC8RO2(i)=fC8RO2(i)+1; 

RxnToReplace = 'C810O = CH3COCH3 + C514O2';
kToReplace = KDEC.*0.84;
ReplaceRxn
i=i+1;
Rnames{i} = 'C810O = C8RO2'; % 1,5-H shift (aldehydic H) to produce C8RO2.
k(:,i) = KDEC.*0.16;
Gstr{i,1} = 'C810O'; 
fC810O(i)=fC810O(i)-1; fC8RO2(i)=fC8RO2(i)+1; 

%% Criegee Intermediate reactions
% the new branching ratios are from Claflin et al., 2018 ACS Earth and Space
% Chem.
RxnToReplace = 'APINOOA = C107O2 + OH';
kToReplace = KDEC.*0.22;
ReplaceRxn

RxnToReplace = 'APINOOA = C109O2 + OH';
kToReplace = KDEC.*0.66;
ReplaceRxn

i=i+1;
Rnames{i} = 'APINOOA = APINAOO';
k(:,i) = KDEC.*0.12;
Gstr{i,1} = 'APINOOA'; 
fAPINOOA(i)=fAPINOOA(i)-1; fAPINAOO(i)=fAPINAOO(i)+1;

RxnToReplace = 'APINOOB = APINBOO';
kToReplace = KDEC.*0.18;
ReplaceRxn

RxnToReplace = 'APINOOB = C96O2 + OH + CO';
kToReplace = KDEC.*0.0;
ReplaceRxn

i=i+1;
Rnames{i} = 'APINOOB = PINONIC';
k(:,i) = KDEC.*0.16;
Gstr{i,1} = 'APINOOB'; 
fAPINOOB(i)=fAPINOOB(i)-1; fPINONIC(i)=fPINONIC(i)+1; 

% the branching ratios of 0.3 vs. 0.7 in the next two rxns are the branching
% ratios for ring-retained vs. ring-opening. See Iyer et al., 2020 Nature
% Comm.
i=i+1;
Rnames{i} = ' APINOOB = PINAL4RO2 + OH ';
k(:,i) = KDEC.*0.66.*(1-fAPINOOB_RB);
Gstr{i,1} = 'APINOOB';
fAPINOOB(i)=fAPINOOB(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)+1; fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = ' APINOOB = APC10H15O4RB + OH ';
k(:,i) = KDEC.*0.66.*fAPINOOB_RB;
Gstr{i,1} = 'APINOOB';
fAPINOOB(i)=fAPINOOB(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)+1; fOH(i)=fOH(i)+1;

% bimolecular rxns for APINAOO. APINBOO rxns are already included in the
% base MCM331
i=i+1;
Rnames{i} = 'APINAOO + CO = PINAL + CO2';
k(:,i) = 1.2E-15;
Gstr{i,1} = 'APINAOO'; Gstr{i,2} = 'CO'; 
fAPINAOO(i)=fAPINAOO(i)-1; fCO(i)=fCO(i)-1; fPINAL(i)=fPINAL(i)+1; 
i=i+1;
Rnames{i} = 'APINAOO + NO = PINAL + NO2';
k(:,i) = 1.2E-15;
Gstr{i,1} = 'APINAOO'; Gstr{i,2} = 'NO'; 
fAPINAOO(i)=fAPINAOO(i)-1; fNO(i)=fNO(i)-1; fPINAL(i)=fPINAL(i)+1; fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APINAOO + NO2 = PINAL + NO3';
k(:,i) = 1.00e-15;
Gstr{i,1} = 'APINAOO'; Gstr{i,2} = 'NO2'; 
fAPINAOO(i)=fAPINAOO(i)-1; fNO2(i)=fNO2(i)-1; fPINAL(i)=fPINAL(i)+1; fNO3(i)=fNO3(i)+1; 
i=i+1;
Rnames{i} = 'APINAOO + SO2 = PINAL + SO3';
k(:,i) = 7.00e-14;
Gstr{i,1} = 'APINAOO'; Gstr{i,2} = 'SO2'; 
fAPINAOO(i)=fAPINAOO(i)-1; fSO2(i)=fSO2(i)-1; fPINAL(i)=fPINAL(i)+1; fSO3(i)=fSO3(i)+1; 
i=i+1;
Rnames{i} = 'APINAOO = PINAL + H2O2';
k(:,i) = 1.40e-17.*H2O;
Gstr{i,1} = 'APINAOO'; 
fAPINAOO(i)=fAPINAOO(i)-1; fPINAL(i)=fPINAL(i)+1; fH2O2(i)=fH2O2(i)+1; 
i=i+1;
Rnames{i} = 'APINAOO = PINONIC';
k(:,i) = 2.00e-18.*H2O;
Gstr{i,1} = 'APINAOO'; 
fAPINAOO(i)=fAPINAOO(i)-1; fPINONIC(i)=fPINONIC(i)+1; 

% Decomposition of SCI to the corresponding RO2. Rate constant based on
% Vereecken et al., PCCP, 2017, 19, 31599.
i=i+1;
Rnames{i} = 'APINAOO = C107O2 + OH';
k(:,i) = 155.*0.75;
Gstr{i,1} = 'APINAOO'; 
fAPINAOO(i)=fAPINAOO(i)-1; fC107O2(i)=fC107O2(i)+1; fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'APINAOO = C109O2 + OH';
k(:,i) = 155.*0.25;
Gstr{i,1} = 'APINAOO'; 
fAPINAOO(i)=fAPINAOO(i)-1; fC109O2(i)=fC109O2(i)+1; fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = ' APINBOO = PINAL4RO2 + OH ';
k(:,i) = 600.*(1-fAPINOOB_RB);
Gstr{i,1} = 'APINBOO';
fAPINBOO(i)=fAPINBOO(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)+1; fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = ' APINBOO = APC10H15O4RB + OH ';
k(:,i) = 600.*fAPINOOB_RB;
Gstr{i,1} = 'APINBOO'; 
fAPINBOO(i)=fAPINBOO(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)+1; fOH(i)=fOH(i)+1;

% SCI + small acids/alcohols/ketones.
i=i+1;
Rnames{i} = 'APINAOO + CH3CO2H = AAHP';
k(:,i) = 5E-10;
Gstr{i,1} = 'APINAOO'; Gstr{i,2} = 'CH3CO2H'; 
fAPINAOO(i)=fAPINAOO(i)-1; fCH3CO2H(i)=fCH3CO2H(i)-1; fAAHP(i)=fAAHP(i)+1;
i=i+1;
Rnames{i} = 'APINAOO + HCOOH = AAHP';
k(:,i) = 5E-10;
Gstr{i,1} = 'APINAOO'; Gstr{i,2} = 'HCOOH'; 
fAPINAOO(i)=fAPINAOO(i)-1; fHCOOH(i)=fHCOOH(i)-1; fAAHP(i)=fAAHP(i)+1;
i=i+1;
Rnames{i} = 'APINAOO + CH3OH = AAHP';
k(:,i) = 4E-14;
Gstr{i,1} = 'APINAOO'; Gstr{i,2} = 'CH3OH'; 
fAPINAOO(i)=fAPINAOO(i)-1; fCH3OH(i)=fCH3OH(i)-1; fAAHP(i)=fAAHP(i)+1;
i=i+1;
Rnames{i} = 'APINAOO + IPROPOL = AAHP';
k(:,i) = 4E-14;
Gstr{i,1} = 'APINAOO'; Gstr{i,2} = 'IPROPOL'; 
fAPINAOO(i)=fAPINAOO(i)-1; fIPROPOL(i)=fIPROPOL(i)-1; fAAHP(i)=fAAHP(i)+1;
i=i+1;
Rnames{i} = 'APINAOO + NBUTOL = AAHP';
k(:,i) = 4E-14;
Gstr{i,1} = 'APINAOO'; Gstr{i,2} = 'NBUTOL'; 
fAPINAOO(i)=fAPINAOO(i)-1; fNBUTOL(i)=fNBUTOL(i)-1; fAAHP(i)=fAAHP(i)+1;
i=i+1;
Rnames{i} = 'APINAOO + BUT2OL = AAHP';
k(:,i) = 4E-14;
Gstr{i,1} = 'APINAOO'; Gstr{i,2} = 'BUT2OL'; 
fAPINAOO(i)=fAPINAOO(i)-1; fBUT2OL(i)=fBUT2OL(i)-1; fAAHP(i)=fAAHP(i)+1;
i=i+1;
Rnames{i} = 'APINAOO + CH3COCH3 = SOZ';
k(:,i) = 4.7E-13;
Gstr{i,1} = 'APINAOO'; Gstr{i,2} = 'CH3COCH3'; 
fAPINAOO(i)=fAPINAOO(i)-1; fCH3COCH3(i)=fCH3COCH3(i)-1; fSOZ(i)=fSOZ(i)+1;

i=i+1;
Rnames{i} = 'APINBOO + CH3CO2H = AAHP';
k(:,i) = 5E-10;
Gstr{i,1} = 'APINBOO'; Gstr{i,2} = 'CH3CO2H'; 
fAPINBOO(i)=fAPINBOO(i)-1; fCH3CO2H(i)=fCH3CO2H(i)-1; fAAHP(i)=fAAHP(i)+1;
i=i+1;
Rnames{i} = 'APINBOO + HCOOH = AAHP';
k(:,i) = 5E-10;
Gstr{i,1} = 'APINBOO'; Gstr{i,2} = 'HCOOH'; 
fAPINBOO(i)=fAPINBOO(i)-1; fHCOOH(i)=fHCOOH(i)-1; fAAHP(i)=fAAHP(i)+1;
i=i+1;
Rnames{i} = 'APINBOO + CH3OH = AAHP';
k(:,i) = 4E-14;
Gstr{i,1} = 'APINBOO'; Gstr{i,2} = 'CH3OH'; 
fAPINBOO(i)=fAPINBOO(i)-1; fCH3OH(i)=fCH3OH(i)-1; fAAHP(i)=fAAHP(i)+1;
i=i+1;
Rnames{i} = 'APINBOO + IPROPOL = AAHP';
k(:,i) = 4E-14;
Gstr{i,1} = 'APINBOO'; Gstr{i,2} = 'IPROPOL'; 
fAPINBOO(i)=fAPINBOO(i)-1; fIPROPOL(i)=fIPROPOL(i)-1; fAAHP(i)=fAAHP(i)+1;
i=i+1;
Rnames{i} = 'APINBOO + NBUTOL = AAHP';
k(:,i) = 4E-14;
Gstr{i,1} = 'APINBOO'; Gstr{i,2} = 'NBUTOL'; 
fAPINBOO(i)=fAPINBOO(i)-1; fNBUTOL(i)=fNBUTOL(i)-1; fAAHP(i)=fAAHP(i)+1;
i=i+1;
Rnames{i} = 'APINBOO + BUT2OL = AAHP';
k(:,i) = 4E-14;
Gstr{i,1} = 'APINBOO'; Gstr{i,2} = 'BUT2OL'; 
fAPINBOO(i)=fAPINBOO(i)-1; fBUT2OL(i)=fBUT2OL(i)-1; fAAHP(i)=fAAHP(i)+1;
i=i+1;
Rnames{i} = 'APINBOO + CH3COCH3 = SOZ';
k(:,i) = 4.7E-13;
Gstr{i,1} = 'APINBOO'; Gstr{i,2} = 'CH3COCH3'; 
fAPINBOO(i)=fAPINBOO(i)-1; fCH3COCH3(i)=fCH3COCH3(i)-1; fSOZ(i)=fSOZ(i)+1;

%% C10H15O4-RO2 autoxidation
i=i+1;
Rnames{i} = 'PINAL4RO2 + autox = PINALPA4RO2';
k(:,i) = kautox_PINAL4RO2; % autoxidation producing C10H15O6. PINAL4P3RO2 is rapidly scrambled to this new RO2.
Gstr{i,1} = 'PINAL4RO2';
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O4RB + autox = APC10H15O6RO2';
k(:,i) = kautox_APC10H15O4RB; % autoxidation producing C10H15O6, rate constant based on Iyer et al., 2020 Nature Comm.
Gstr{i,1} = 'APC10H15O4RB';
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)+1;
i=i+1;
Rnames{i} = 'C109O2 + autox = PINAL10P1RO2';
k(:,i) = kautox_C109O2; % autoxidation producing C10H15O6.
Gstr{i,1} = 'C109O2';
fC109O2(i)=fC109O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)+1;
i=i+1;
Rnames{i} = 'C107O2 + autox = PINALPA1RO2';
k(:,i) = kautox_C107O2; % autoxidation producing C10H15O6. PINAL1P3RO2 is rapidly scrambled to this new RO2.
Gstr{i,1} = 'C107O2';
fC107O2(i)=fC107O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)+1;
%% C10H15O4-RO2 bimolecular rxns
% need to consider revise C107O2 and C109O2 + HO2 rxns.
i=i+1;
Rnames{i} = 'PINAL4RO2 + HO2 = PINAL4P'; % C10H16O4
k(:,i) = KRO2HO2.*0.914.*fROOH_sec;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'HO2';
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fHO2(i)=fHO2(i)-1;fPINAL4P(i)=fPINAL4P(i)+1;
i=i+1;
Rnames{i} = 'PINAL4RO2 + HO2 = PINAL4RO + OH';
k(:,i) = KRO2HO2.*0.914.*(1.0-fROOH_sec);
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'HO2';
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fHO2(i)=fHO2(i)-1;fPINAL4RO(i)=fPINAL4RO(i)+1; fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O4RB + HO2 = APC10H16O4'; % C10H16O4
k(:,i) = KRO2HO2.*0.914.*fROOH_tert;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'HO2';
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H16O4(i)=fAPC10H16O4(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O4RB + HO2 = APC10H15O3RO + OH';
k(:,i) = KRO2HO2.*0.914.*(1.0-fROOH_tert);
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'HO2';
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H15O3RO(i)=fAPC10H15O3RO(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'PINAL4RO2 + NO = PINAL4RO + NO2';
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'NO';
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fNO(i)=fNO(i)-1;fPINAL4RO(i)=fPINAL4RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'PINAL4RO2 + NO = APC10H15NO5';
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'NO';
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15NO5(i)=fAPC10H15NO5(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O4RB + NO = APC10H15O3RO + NO2';
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'NO';
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fNO(i)=fNO(i)-1;fAPC10H15O3RO(i)=fAPC10H15O3RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O4RB + NO = APC10H15NO5';
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'NO';
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fNO(i)=fNO(i)-1;fAPC10H15NO5(i)=fAPC10H15NO5(i)+1;
i=i+1;
Rnames{i} = 'PINAL4RO2 + NO3 = PINAL4RO + NO2';
k(:,i) = KRO2NO3;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'NO3';
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fNO3(i)=fNO3(i)-1;fPINAL4RO(i)=fPINAL4RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O4RB + NO3 = APC10H15O3RO + NO2';
k(:,i) = KRO2NO3;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'NO';
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fNO3(i)=fNO3(i)-1;fAPC10H15O3RO(i)=fAPC10H15O3RO(i)+1;fNO2(i)=fNO2(i)+1;

i=i+1;
Rnames{i} = 'PINAL4RO2 + RO2 = PINAL4H'; % C10H16O3
k(:,i) = KRO2RO2_fast.*bROH_ps;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'RO2';
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fPINAL4H(i)=fPINAL4H(i)+1;
i=i+1;
Rnames{i} = 'PINAL4RO2 + RO2 = PINAL4C'; % C10H14O3
k(:,i) = KRO2RO2_fast.*bROH_ps;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'RO2';
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fPINAL4C(i)=fPINAL4C(i)+1;
i=i+1;
Rnames{i} = 'PINAL4RO2 + RO2 = PINAL4RO'; % C10H15O3-RO
k(:,i) = KRO2RO2_fast.*bRO_ps;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'RO2';
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fPINAL4RO(i)=fPINAL4RO(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O4RB + RO2 = APC10H16O3'; % C10H16O3
k(:,i) = KRO2RO2_fast.*bROH_ta;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'RO2';
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fAPC10H16O3(i)=fAPC10H16O3(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O4RB + RO2 = APC10H15O3RO'; % C10H15O3-RO
k(:,i) = KRO2RO2_fast.*bRO_ta;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'RO2';
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fAPC10H15O3RO(i)=fAPC10H15O3RO(i)+1;
i=i+1;
Rnames{i} = 'PINAL4RO + decomp = NORPINAL + HO2 + CO'; % RO decomposition. dominant pathway for this RO based on SAR
k(:,i) = KDEC.*0.91;
Gstr{i,1} = 'PINAL4RO';
fPINAL4RO(i)=fPINAL4RO(i)-1; fNORPINAL(i)=fNORPINAL(i)+1; fHO2(i)=fHO2(i)+1; fCO(i)=fCO(i)+1;
i=i+1;
Rnames{i} = 'PINAL4RO + isom = PINAL4H1RO2'; % RO 1,5-H shift isomerization, producing C10H15O5-RO2, RO2 likely on 1 position, branching ratio vs. decomp based on SAR
k(:,i) = KDEC.*0.09;
Gstr{i,1} = 'PINAL4RO';
fPINAL4RO(i)=fPINAL4RO(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O3RO + decomp = NORPINAL + HO2 + CO'; % RO decomposition, dominant pathway for this RO based on SAR.
k(:,i) = KDEC.*0.91;
Gstr{i,1} = 'APC10H15O3RO';
fAPC10H15O3RO(i)=fAPC10H15O3RO(i)-1; fNORPINAL(i)=fNORPINAL(i)+1; fHO2(i)=fHO2(i)+1; fCO(i)=fCO(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O3RO + isom = APC10H15O5RO2'; % RO isomerization, producing another RO2, C10H15O5
k(:,i) = KDEC.*0.09;
Gstr{i,1} = 'APC10H15O3RO';
fAPC10H15O3RO(i)=fAPC10H15O3RO(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)+1;

%% C10H15O5-RO2 autoxidation
i=i+1;
Rnames{i} = 'PINAL4H1RO2 + autox = APC10H15O7RO2';
k(:,i) = kautox_PINAL4H1RO2; % autoxidation producing C10H15O7, rate constant based on SAR.
Gstr{i,1} = 'PINAL4H1RO2';
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O5RO2 + autox = APC10H15O7RO2';
k(:,i) = kautox_C10H15O3O2; % autoxidation producing C10H15O7, rate constant unknown. Using PRAM value for C10H15O3O2
Gstr{i,1} = 'APC10H15O5RO2';
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)+1;
i=i+1;
Rnames{i} = 'C108O2 + autox = APC10H15O7RO2';
k(:,i) = kautox_C108O2; % autoxidation producing C10H15O7, rate constant based on SAR.
Gstr{i,1} = 'C108O2';
fC108O2(i)=fC108O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)+1;
i=i+1;
Rnames{i} = 'C920CO3 + autox = APC10H15O7RO2';
k(:,i) = kautox_C920CO3; % autoxidation producing C10H15O7, rate constant based on SAR.
Gstr{i,1} = 'C920CO3';
fC920CO3(i)=fC920CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)+1;

%% C10H15O5-RO2 bimolecular rxns
i=i+1;
Rnames{i} = 'PINAL4H1RO2 + HO2 = PINAL4H1P'; % C10H16O5
k(:,i) = KRO2HO2.*0.914.*fROOH_tert;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'HO2';
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fHO2(i)=fHO2(i)-1;fPINAL4H1P(i)=fPINAL4H1P(i)+1;
i=i+1;
Rnames{i} = 'PINAL4H1RO2 + HO2 = PINAL4H1RO + OH'; % C10H15O4-RO
k(:,i) = KRO2HO2.*0.914.*(1.0-fROOH_tert);
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'HO2';
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fHO2(i)=fHO2(i)-1;fPINAL4H1RO(i)=fPINAL4H1RO(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O5RO2 + HO2 = APC10H16O5'; % C10H16O5. This RO2 is assumed to be 50% secondary and 50% tertiary
k(:,i) = KRO2HO2.*0.914.*(0.5*fROOH_sec+0.5.*fROOH_tert);
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'HO2';
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H16O5(i)=fAPC10H16O5(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O5RO2 + HO2 = APC10H15O4RO + OH'; % C10H15O4-RO
k(:,i) = KRO2HO2.*0.914.*(0.5*(1.0-fROOH_sec)+0.5.*(1.0-fROOH_tert));
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'HO2';
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H15O4RO(i)=fAPC10H15O4RO(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'PINAL4H1RO2 + NO = PINAL4H1RO + NO2';
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'NO';
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fNO(i)=fNO(i)-1;fPINAL4H1RO(i)=fPINAL4H1RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'PINAL4H1RO2 + NO = APC10H15NO6';
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'NO';
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15NO6(i)=fAPC10H15NO6(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O5RO2 + NO = APC10H15O4RO + NO2';
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'NO';
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15O4RO(i)=fAPC10H15O4RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O5RO2 + NO = APC10H15NO6';
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'NO';
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15NO6(i)=fAPC10H15NO6(i)+1;

i=i+1;
Rnames{i} = 'PINAL4H1RO2 + NO3 = PINAL4H1RO + NO2';
k(:,i) = KRO2NO3;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'NO3';
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fNO3(i)=fNO3(i)-1;fPINAL4H1RO(i)=fPINAL4H1RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O5RO2 + NO3 = APC10H15O4RO + NO2';
k(:,i) = KRO2NO3;
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'NO3';
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fNO3(i)=fNO3(i)-1;fAPC10H15O4RO(i)=fAPC10H15O4RO(i)+1;fNO2(i)=fNO2(i)+1;

i=i+1;
Rnames{i} = 'PINAL4H1RO2 + RO2 = PINAL1H4H'; % C10H16O4
k(:,i) = KRO2RO2_fast.*bROH_ta;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'RO2';
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fPINAL1H4H(i)=fPINAL1H4H(i)+1;
i=i+1;
Rnames{i} = 'PINAL4H1RO2 + RO2 = PINAL4H1RO'; % C10H15O4-RO
k(:,i) = KRO2RO2_fast.*bRO_ta;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'RO2';
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fPINAL4H1RO(i)=fPINAL4H1RO(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O5RO2 + RO2 = APC10H16O4'; % C10H16O4. For lumped HOM-RO2, 50% secondary and 50% tertiary
k(:,i) = KRO2RO2_fast.*(0.5.*bROH_ps + 0.5.*bROH_ta);
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fAPC10H16O4(i)=fAPC10H16O4(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O5RO2 + RO2 = APC10H14O4'; % C10H14O4. For lumped HOM-RO2, 50% secondary and 50% tertiary
k(:,i) = KRO2RO2_fast.*0.5.*bROH_ps;
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fAPC10H14O4(i)=fAPC10H14O4(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O5RO2 + RO2 = APC10H15O4RO'; % For lumped HOM-RO2, 50% secondary and 50% tertiary
k(:,i) = KRO2RO2_fast.*(0.5.*bRO_ps + 0.5.*bRO_ta);
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fAPC10H15O4RO(i)=fAPC10H15O4RO(i)+1;

i=i+1;
Rnames{i} = 'PINAL4H1RO + decomp = APC10H15O6RO2'; % RO decomposition, ring-opening, similar to C107O, C10H15O6-RO2. dominant pathway for this RO based on SAR
k(:,i) = KDEC;
Gstr{i,1} = 'PINAL4H1RO';
fPINAL4H1RO(i)=fPINAL4H1RO(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O4RO + decomp = APC10H15O6RO2'; % RO decomposition, ring-opening. dominant pathway for this RO based on SAR
k(:,i) = KDEC;
Gstr{i,1} = 'APC10H15O4RO';
fAPC10H15O4RO(i)=fAPC10H15O4RO(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)+1;

%% C10H15O6-RO2 autoxidation
i=i+1;
Rnames{i} = 'PINALPA4RO2 + autox = APC10H15O8RO2';
k(:,i) = kautox_PINALPA4RO2; % autoxidation producing C10H15O8, rate constant based on SAR
Gstr{i,1} = 'PINALPA4RO2';
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)+1;
i=i+1;
Rnames{i} = 'PINALPA1RO2 + autox = APC10H15O8RO2';
k(:,i) = kautox_PINALPA1RO2; % autoxidation producing C10H15O8, rate constant based on SAR
Gstr{i,1} = 'PINALPA1RO2';
fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)+1;
i=i+1;
Rnames{i} = 'PINAL10P1RO2 + autox = APC10H15O8RO2';
k(:,i) = kautox_PINAL10P1RO2; % autoxidation producing C10H15O8, rate constant based on SAR
Gstr{i,1} = 'PINAL10P1RO2';
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O6RO2 + autox = APC10H15O8RO2';
k(:,i) = kautox_APC10H15O6RO2; % autoxidation producing C10H15O8, rate constant based on SAR
Gstr{i,1} = 'APC10H15O6RO2';
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)+1;

i=i+1;
Rnames{i} = 'PINAL10P1RO2 + term = APC10H14O5 + OH';
k(:,i) = 0.004; % RO2 intramolecular termination, C10H14O5. Not important based on SAR
Gstr{i,1} = 'PINAL10P1RO2';
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fAPC10H14O5(i)=fAPC10H14O5(i)+1;fOH(i)=fOH(i)+1;

%% C10H15O6-RO2 bimolecular rxns
i=i+1;
Rnames{i} = 'PINALPA4RO2 + HO2 = APC10H16O6'; % secRO2 + HO2, C10H16O6
k(:,i) = KRO2HO2.*0.914.*fROOH_sec;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'HO2';
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H16O6(i)=fAPC10H16O6(i)+1;
i=i+1;
Rnames{i} = 'PINALPA4RO2 + HO2 = PINALPA4RO + OH'; % secRO2 + HO2 = RO + OH
k(:,i) = KRO2HO2.*0.914.*(1.0-fROOH_sec);
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'HO2';
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fHO2(i)=fHO2(i)-1;fPINALPA4RO(i)=fPINALPA4RO(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'PINALPA1RO2 + HO2 = APC10H16O6'; % tertRO2 + HO2, C10H16O6
k(:,i) = KRO2HO2.*0.914.*fROOH_tert;
Gstr{i,1} = 'PINALPA1RO2'; Gstr{i,2} = 'HO2';
fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H16O6(i)=fAPC10H16O6(i)+1;
i=i+1;
Rnames{i} = 'PINALPA1RO2 + HO2 = PINALPA1RO + OH'; % tertRO2 + HO2 = RO + OH
k(:,i) = KRO2HO2.*0.914.*(1.0-fROOH_tert);
Gstr{i,1} = 'PINALPA1RO2'; Gstr{i,2} = 'HO2';
fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fHO2(i)=fHO2(i)-1;fPINALPA1RO(i)=fPINALPA1RO(i)+1;fOH(i)=fOH(i)+1;

i=i+1;
Rnames{i} = 'PINAL10P1RO2 + HO2 = APC10H16O6'; % RO2 + HO2, C10H16O6
k(:,i) = KRO2HO2.*0.914.*fROOH_tert;
Gstr{i,1} = 'PINAL10P1RO2'; Gstr{i,2} = 'HO2';
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H16O6(i)=fAPC10H16O6(i)+1;
i=i+1;
Rnames{i} = 'PINAL10P1RO2 + HO2 = PINAL10P1RO + OH'; % RO2 + HO2 = RO + OH
k(:,i) = KRO2HO2.*0.914.*(1.0-fROOH_tert);
Gstr{i,1} = 'PINAL10P1RO2'; Gstr{i,2} = 'HO2';
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fHO2(i)=fHO2(i)-1;fPINAL10P1RO(i)=fPINAL10P1RO(i)+1;fOH(i)=fOH(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O6RO2 + HO2 = APC10H16O6'; % RO2 + HO2, APC10H16O6
k(:,i) = KRO2HO2.*0.914.*fROOH_sec;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'HO2';
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H16O6(i)=fAPC10H16O6(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O6RO2 + HO2 = APC10H15O5RO + OH'; % RO2 + HO2 = RO + OH
k(:,i) = KRO2HO2.*0.914.*(1.0-fROOH_sec);
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'HO2';
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H15O5RO(i)=fAPC10H15O5RO(i)+1;fOH(i)=fOH(i)+1;

i=i+1;
Rnames{i} = 'PINALPA4RO2 + NO = PINALPA4RO + NO2';
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'NO';
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fNO(i)=fNO(i)-1;fPINALPA4RO(i)=fPINALPA4RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'PINALPA4RO2 + NO = APC10H15NO7';
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'NO';
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15NO7(i)=fAPC10H15NO7(i)+1;
i=i+1;
Rnames{i} = 'PINALPA1RO2 + NO = PINALPA1RO + NO2';
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'PINALPA1RO2'; Gstr{i,2} = 'NO';
fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fNO(i)=fNO(i)-1;fPINALPA1RO(i)=fPINALPA1RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'PINALPA1RO2 + NO = APC10H15NO7';
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'PINALPA1RO2'; Gstr{i,2} = 'NO';
fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15NO7(i)=fAPC10H15NO7(i)+1;
i=i+1;
Rnames{i} = 'PINAL10P1RO2 + NO = PINAL10P1RO + NO2';
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'PINAL10P1RO2'; Gstr{i,2} = 'NO';
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fNO(i)=fNO(i)-1;fPINAL10P1RO(i)=fPINAL10P1RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'PINAL10P1RO2 + NO = APC10H15NO7';
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'PINAL10P1RO2'; Gstr{i,2} = 'NO';
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15NO7(i)=fAPC10H15NO7(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O6RO2 + NO = APC10H15O5RO + NO2';
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'NO';
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15O5RO(i)=fAPC10H15O5RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O6RO2 + NO = APC10H15NO7';
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'NO';
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15NO7(i)=fAPC10H15NO7(i)+1;

i=i+1;
Rnames{i} = 'PINALPA4RO2 + NO3 = PINALPA4RO + NO2';
k(:,i) = KRO2NO3;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'NO3';
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fNO3(i)=fNO3(i)-1;fPINALPA4RO(i)=fPINALPA4RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'PINALPA1RO2 + NO3 = PINALPA1RO + NO2';
k(:,i) = KRO2NO3;
Gstr{i,1} = 'PINALPA1RO2'; Gstr{i,2} = 'NO3';
fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fNO3(i)=fNO3(i)-1;fPINALPA1RO(i)=fPINALPA1RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'PINAL10P1RO2 + NO3 = PINAL10P1RO + NO2';
k(:,i) = KRO2NO3;
Gstr{i,1} = 'PINAL10P1RO2'; Gstr{i,2} = 'NO3';
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fNO3(i)=fNO3(i)-1;fPINAL10P1RO(i)=fPINAL10P1RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O6RO2 + NO3 = APC10H15O5RO + NO2';
k(:,i) = KRO2NO3;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'NO3';
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fNO3(i)=fNO3(i)-1;fAPC10H15O5RO(i)=fAPC10H15O5RO(i)+1;fNO2(i)=fNO2(i)+1;

i=i+1;
Rnames{i} = 'PINALPA4RO2 + RO2 = APC10H16O5'; % C10H16O5
k(:,i) = KRO2RO2_fast.*bROH_ps;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'RO2';
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fAPC10H16O5(i)=fAPC10H16O5(i)+1;
i=i+1;
Rnames{i} = 'PINALPA4RO2 + RO2 = APC10H14O5'; % C10H14O5
k(:,i) = KRO2RO2_fast.*bROH_ps;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'RO2';
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fAPC10H14O5(i)=fAPC10H14O5(i)+1;
i=i+1;
Rnames{i} = 'PINALPA4RO2 + RO2 = PINALPA4RO'; % C10H15O5-RO
k(:,i) = KRO2RO2_fast.*bRO_ps;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'RO2';
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fPINALPA4RO(i)=fPINALPA4RO(i)+1;
i=i+1;
Rnames{i} = 'PINALPA1RO2 + RO2 = APC10H16O5'; % C10H16O5
k(:,i) = KRO2RO2_fast.*bROH_ta;
Gstr{i,1} = 'PINALPA1RO2'; Gstr{i,2} = 'RO2';
fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fAPC10H16O5(i)=fAPC10H16O5(i)+1;
i=i+1;
Rnames{i} = 'PINALPA1RO2 + RO2 = PINALPA1RO'; % C10H15O5-RO
k(:,i) = KRO2RO2_fast.*bRO_ta;
Gstr{i,1} = 'PINALPA1RO2'; Gstr{i,2} = 'RO2';
fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fPINALPA1RO(i)=fPINALPA1RO(i)+1;
i=i+1;
Rnames{i} = 'PINAL10P1RO2 + RO2 = APC10H16O5'; % C10H16O5
k(:,i) = KRO2RO2_fast.*bROH_ta;
Gstr{i,1} = 'PINAL10P1RO2'; Gstr{i,2} = 'RO2';
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fAPC10H16O5(i)=fAPC10H16O5(i)+1;
i=i+1;
Rnames{i} = 'PINAL10P1RO2 + RO2 = PINAL10P1RO'; % C10H15O5-RO
k(:,i) = KRO2RO2_fast.*bRO_ta;
Gstr{i,1} = 'PINAL10P1RO2'; Gstr{i,2} = 'RO2';
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fPINAL10P1RO(i)=fPINAL10P1RO(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O6RO2 + RO2 = APC10H16O5'; % C10H16O5. secondary RO2
k(:,i) = KRO2RO2_fast.*bROH_ps;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fAPC10H16O5(i)=fAPC10H16O5(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O6RO2 + RO2 = APC10H14O5'; % C10H14O5
k(:,i) = KRO2RO2_fast.*bROH_ps;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fAPC10H14O5(i)=fAPC10H14O5(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O6RO2 + RO2 = APC10H15O5RO'; % C10H15O5-RO
k(:,i) = KRO2RO2_fast.*bRO_ps;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fAPC10H15O5RO(i)=fAPC10H15O5RO(i)+1;

i=i+1;
Rnames{i} = 'PINALPA4RO + decomp = NORPINAL + OH'; % RO decomposition. Dominant pathway for this RO based on SAR.
k(:,i) = KDEC;
Gstr{i,1} = 'PINALPA4RO';
fPINALPA4RO(i)=fPINALPA4RO(i)-1; fNORPINAL(i)=fNORPINAL(i)+1; fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'PINALPA1RO + decomp = APC10H15O7RO2'; % RO ring-opening. Dominant pathway for this RO based on SAR.
k(:,i) = KDEC;
Gstr{i,1} = 'PINALPA1RO';
fPINALPA1RO(i)=fPINALPA1RO(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)+1;
i=i+1;
Rnames{i} = 'PINAL10P1RO + decomp = APC10H15O7RO2'; % RO decomposition. producing a ring-opening RO2, APC10H15O7RO2 Dominant pathway for this RO based on SAR
k(:,i) = KDEC;
Gstr{i,1} = 'PINAL10P1RO';
fPINAL10P1RO(i)=fPINAL10P1RO(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O5RO + decomp = C9C'; % endoperoxide RO decomposition, dominant pathway for this RO based on SAR. Producing C9H14O4 (isomer of pinic acid, without exchangable H)
k(:,i) = KDEC.*0.91;
Gstr{i,1} = 'APC10H15O5RO';
fAPC10H15O5RO(i)=fAPC10H15O5RO(i)-1; fC9C(i)=fC9C(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O5RO + isom = APC10H15O7RO2'; % endoperoxide RO 1,5-H isomerization
k(:,i) = KDEC.*0.09;
Gstr{i,1} = 'APC10H15O5RO';
fAPC10H15O5RO(i)=fAPC10H15O5RO(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)+1;

%% C10H15O7-RO2 autoxidation
i=i+1;
Rnames{i} = 'APC10H15O7RO2 + autox = APC10H15O9RO2';
k(:,i) = kautox_APC10H15O7RO2; % autoxidation producing C10H15O9, rate constant based on SAR for two 50:50 isomers.
Gstr{i,1} = 'APC10H15O7RO2';
fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)+1;
%% C10H15O7-RO2 bimolecular rxns
i=i+1;
Rnames{i} = 'APC10H15O7RO2 + HO2 = APC10H16O7'; % RO2 + HO2, C10H16O7. tertiary RO2
k(:,i) = KRO2HO2.*0.914.*fROOH_tert;
Gstr{i,1} = 'APC10H15O7RO2'; Gstr{i,2} = 'HO2';
fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H16O7(i)=fAPC10H16O7(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O7RO2 + HO2 = APC10H15O6RO + OH'; % RO2 + HO2 = RO + OH.  tertiary RO2
k(:,i) = KRO2HO2.*0.914.*(1.0-fROOH_tert);
Gstr{i,1} = 'APC10H15O7RO2'; Gstr{i,2} = 'HO2';
fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H15O6RO(i)=fAPC10H15O6RO(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O7RO2 + NO = APC10H15O6RO + NO2';
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'APC10H15O7RO2'; Gstr{i,2} = 'NO';
fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15O6RO(i)=fAPC10H15O6RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O7RO2 + NO = APC10H15NO8';
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'APC10H15O7RO2'; Gstr{i,2} = 'NO';
fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15NO8(i)=fAPC10H15NO8(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O7RO2 + NO3 = APC10H15O6RO + NO2';
k(:,i) = KRO2NO3;
Gstr{i,1} = 'APC10H15O7RO2'; Gstr{i,2} = 'NO3';
fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fNO3(i)=fNO3(i)-1;fAPC10H15O6RO(i)=fAPC10H15O6RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O7RO2 + RO2 = APC10H16O6'; % C10H16O6
k(:,i) = KRO2RO2_fast.*bROH_ta;
Gstr{i,1} = 'APC10H15O7RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fAPC10H16O6(i)=fAPC10H16O6(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O7RO2 + RO2 = APC10H15O6RO'; % C10H15O6-RO
k(:,i) = KRO2RO2_fast.*bRO_ta;
Gstr{i,1} = 'APC10H15O7RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fAPC10H15O6RO(i)=fAPC10H15O6RO(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O6RO + decomp = C7RO2 + CH3COCH3'; % RO decomposition. The branching ratio between the two pathways are guessed.
k(:,i) = KDEC.*0.9;
Gstr{i,1} = 'APC10H15O6RO';
fAPC10H15O6RO(i)=fAPC10H15O6RO(i)-1; fC7RO2(i)=fC7RO2(i)+1;fCH3COCH3(i)=fCH3COCH3(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O6RO + isom = APC10H15O8RO2'; % tert RO 1,5-H aldehydic isomerization, producing another RO2, C10H15O8
k(:,i) = KDEC.*0.1;
Gstr{i,1} = 'APC10H15O6RO';
fAPC10H15O6RO(i)=fAPC10H15O6RO(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)+1;

%% C10H15O8-RO2 autoxidation
i=i+1;
Rnames{i} = 'APC10H15O8RO2 + autox = APC10H15O10RO2';
k(:,i) = kautox_APC10H15O8RO2; % autoxidation producing C10H15O10, rate constant by SAR.
Gstr{i,1} = 'APC10H15O8RO2';
fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)+1;

%% C10H15O8-RO2 bimolecular rxns
i=i+1;
Rnames{i} = 'APC10H15O8RO2 + HO2 = APC10H16O8'; % RO2 + HO2, C10H16O8
k(:,i) = KRO2HO2.*0.914.*(fsec.*fROOH_sec+ftert.*fROOH_tert);
Gstr{i,1} = 'APC10H15O8RO2'; Gstr{i,2} = 'HO2';
fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H16O8(i)=fAPC10H16O8(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O8RO2 + HO2 = APC10H15O7RO + OH'; % RO2 + HO2 = RO + OH
k(:,i) = KRO2HO2.*0.914.*(ftert.*(1.0-fROOH_tert)+fsec.*(1.0-fROOH_sec));
Gstr{i,1} = 'APC10H15O8RO2'; Gstr{i,2} = 'HO2';
fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H15O7RO(i)=fAPC10H15O7RO(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O8RO2 + NO = APC10H15O7RO + NO2';
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'APC10H15O8RO2'; Gstr{i,2} = 'NO';
fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15O7RO(i)=fAPC10H15O7RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O8RO2 + NO = APC10H15NO9';
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'APC10H15O8RO2'; Gstr{i,2} = 'NO';
fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15NO9(i)=fAPC10H15NO9(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O8RO2 + NO3 = APC10H15O7RO + NO2';
k(:,i) = KRO2NO3;
Gstr{i,1} = 'APC10H15O8RO2'; Gstr{i,2} = 'NO3';
fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fNO3(i)=fNO3(i)-1;fAPC10H15O7RO(i)=fAPC10H15O7RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O8RO2 + RO2 = APC10H16O7'; % C10H16O7
k(:,i) = KRO2RO2_fast.*(ftert.*bROH_ta+fsec.*bROH_ps); 
Gstr{i,1} = 'APC10H15O8RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fAPC10H16O7(i)=fAPC10H16O7(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O8RO2 + RO2 = APC10H14O7'; % C10H14O7
k(:,i) = KRO2RO2_fast.*fsec.*bROH_ps; 
Gstr{i,1} = 'APC10H15O8RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fAPC10H14O7(i)=fAPC10H14O7(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O8RO2 + RO2 = APC10H15O7RO'; % C10H15O7-RO
k(:,i) = KRO2RO2_fast.*(ftert.*bRO_ta+fsec.*bRO_ps);
Gstr{i,1} = 'APC10H15O8RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fAPC10H15O7RO(i)=fAPC10H15O7RO(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O7RO + decomp = C9OOH + HO2'; % RO decomposition
k(:,i) = KDEC.*0.7;
Gstr{i,1} = 'APC10H15O7RO';
fAPC10H15O7RO(i)=fAPC10H15O7RO(i)-1; fC9OOH(i)=fC9OOH(i)+1;fHO2(i)=fHO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O7RO + decomp = C9OOH + OH'; % RO decomposition
k(:,i) = KDEC.*0.3;
Gstr{i,1} = 'APC10H15O7RO';
fAPC10H15O7RO(i)=fAPC10H15O7RO(i)-1; fC9OOH(i)=fC9OOH(i)+1; fOH(i)=fOH(i)+1;

%% C10H15O9-RO2 autoxidation
i=i+1;
Rnames{i} = 'APC10H15O9RO2 + term = APC10H14O8 + OH';
k(:,i) = 0.0; % RO2 intramolecular termination, C10H14O8. Not important based on SAR.
Gstr{i,1} = 'APC10H15O9RO2';
fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fAPC10H14O8(i)=fAPC10H14O8(i)+1;fOH(i)=fOH(i)+1;
%% C10H15O9-RO2 bimolecular rxns
i=i+1;
Rnames{i} = 'APC10H15O9RO2 + HO2 = APC10H16O9'; % RO2 + HO2, C10H16O9. Assuming 30% of this RO2 is tertiary RO2.
k(:,i) = KRO2HO2.*0.914.*(0.3.*fROOH_tert+0.7.*fROOH_sec);
Gstr{i,1} = 'APC10H15O9RO2'; Gstr{i,2} = 'HO2';
fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H16O9(i)=fAPC10H16O9(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O9RO2 + HO2 = APC10H15O8RO + OH'; % RO2 + HO2 = RO + OH
k(:,i) = KRO2HO2.*0.914.*(0.3.*(1.0-fROOH_tert)+0.7.*(1.0-fROOH_sec));
Gstr{i,1} = 'APC10H15O9RO2'; Gstr{i,2} = 'HO2';
fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H15O8RO(i)=fAPC10H15O8RO(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O9RO2 + NO = APC10H15O8RO + NO2';
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'APC10H15O9RO2'; Gstr{i,2} = 'NO';
fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15O8RO(i)=fAPC10H15O8RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O9RO2 + NO = APC10H15NO10';
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'APC10H15O9RO2'; Gstr{i,2} = 'NO';
fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15NO10(i)=fAPC10H15NO10(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O9RO2 + NO3 = APC10H15O8RO + NO2';
k(:,i) = KRO2NO3;
Gstr{i,1} = 'APC10H15O9RO2'; Gstr{i,2} = 'NO3';
fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fNO3(i)=fNO3(i)-1;fAPC10H15O8RO(i)=fAPC10H15O8RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O9RO2 + RO2 = APC10H16O8'; % C10H16O8
k(:,i) = KRO2RO2_fast.*(0.3.*bROH_ta+0.7.*bROH_ps);
Gstr{i,1} = 'APC10H15O9RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fAPC10H16O8(i)=fAPC10H16O8(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O9RO2 + RO2 = APC10H14O8'; % C10H14O8
k(:,i) = KRO2RO2_fast.*(0.7.*bROH_ps);
Gstr{i,1} = 'APC10H15O9RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fAPC10H14O8(i)=fAPC10H14O8(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O9RO2 + RO2 = APC10H15O8RO'; % C10H15O8-RO
k(:,i) = KRO2RO2_fast.*(0.3.*bRO_ta+0.7.*bRO_ps);
Gstr{i,1} = 'APC10H15O9RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fAPC10H15O8RO(i)=fAPC10H15O8RO(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O8RO + decomp = C7RO2 + CH3COCH3'; % RO decomposition. Branching ratios are guessed.
k(:,i) = KDEC.*0.3;
Gstr{i,1} = 'APC10H15O8RO';
fAPC10H15O8RO(i)=fAPC10H15O8RO(i)-1; fC7RO2(i)=fC7RO2(i)+1; fCH3COCH3(i)=fCH3COCH3(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O8RO + decomp = C7OOH + CH3COCO3'; % RO decomposition
k(:,i) = KDEC.*0.7;
Gstr{i,1} = 'APC10H15O8RO';
fAPC10H15O8RO(i)=fAPC10H15O8RO(i)-1; fC7OOH(i)=fC7OOH(i)+1; fCH3COCO3(i)=fCH3COCO3(i)+1;

%% C10H15O10-RO2 autoxidation
i=i+1;
Rnames{i} = 'APC10H15O10RO2 + autox = APC10H14O9 + OH';
k(:,i) = kautox_APC10H15O10RO2; % RO2 intramolecular termination, C10H14O9
Gstr{i,1} = 'APC10H15O10RO2';
fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fAPC10H14O9(i)=fAPC10H14O9(i)+1;fOH(i)=fOH(i)+1;

%% C10H15O10-RO2 bimolecular rxns
i=i+1;
Rnames{i} = 'APC10H15O10RO2 + HO2 = APC10H16O10'; % RO2 + HO2, C10H16O10
k(:,i) = KRO2HO2.*0.914.*fROOH_sec;
Gstr{i,1} = 'APC10H15O10RO2'; Gstr{i,2} = 'HO2';
fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H16O10(i)=fAPC10H16O10(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O10RO2 + HO2 = APC10H15O9RO + OH'; % RO2 + HO2 = RO + OH
k(:,i) = KRO2HO2.*0.914.*(1.0-fROOH_sec);
Gstr{i,1} = 'APC10H15O10RO2'; Gstr{i,2} = 'HO2';
fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fHO2(i)=fHO2(i)-1;fAPC10H15O9RO(i)=fAPC10H15O9RO(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O10RO2 + NO = APC10H15O9RO + NO2';
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'APC10H15O10RO2'; Gstr{i,2} = 'NO';
fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15O9RO(i)=fAPC10H15O9RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O10RO2 + NO = APC10H15NO11';
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'APC10H15O10RO2'; Gstr{i,2} = 'NO';
fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fNO(i)=fNO(i)-1;fAPC10H15NO11(i)=fAPC10H15NO11(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O10RO2 + NO3 = APC10H15O9RO + NO2';
k(:,i) = KRO2NO3;
Gstr{i,1} = 'APC10H15O10RO2'; Gstr{i,2} = 'NO3';
fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fNO3(i)=fNO3(i)-1;fAPC10H15O9RO(i)=fAPC10H15O9RO(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O10RO2 + RO2 = APC10H16O9'; % C10H16O9
k(:,i) = KRO2RO2_fast.*bROH_ps;
Gstr{i,1} = 'APC10H15O10RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fAPC10H16O9(i)=fAPC10H16O9(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O10RO2 + RO2 = APC10H14O9'; % C10H14O9
k(:,i) = KRO2RO2_fast.*bROH_ps;
Gstr{i,1} = 'APC10H15O10RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fAPC10H14O9(i)=fAPC10H14O9(i)+1;
i=i+1;
Rnames{i} = 'APC10H15O10RO2 + RO2 = APC10H15O9RO'; % C10H15O9-RO
k(:,i) = KRO2RO2_fast.*bRO_ps;
Gstr{i,1} = 'APC10H15O10RO2'; Gstr{i,2} = 'RO2';
fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fAPC10H15O9RO(i)=fAPC10H15O9RO(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O9RO + decomp = C9OOH + OH'; % RO decomposition
k(:,i) = KDEC;
Gstr{i,1} = 'APC10H15O9RO';
fAPC10H15O9RO(i)=fAPC10H15O9RO(i)-1; fC9OOH(i)=fC9OOH(i)+1; fOH(i)=fOH(i)+1;

%% Lumped fragmentation RO2 rxns
i=i+1;
Rnames{i} = 'C9RO2 + autox = C9RO2';
k(:,i) = kautox_gen; % autoxidation (more oxidized C9RO2)
Gstr{i,1} = 'C9RO2'; 
fC9RO2(i)=fC9RO2(i)-1; fC9RO2(i)=fC9RO2(i)+1;
i=i+1;
Rnames{i} = 'C9RO2 + HO2 = C9OOH'; % RO2 + HO2 = ROOH
k(:,i) = KRO2HO2.*0.914.*(0.5.*fROOH_sec + 0.5.*fROOH_tert);
Gstr{i,1} = 'C9RO2'; Gstr{i,2} = 'HO2'; 
fC9RO2(i)=fC9RO2(i)-1; fHO2(i)=fHO2(i)-1;fC9OOH(i)=fC9OOH(i)+1;
i=i+1;
Rnames{i} = 'C9RO2 + HO2 = fragRO2 + OH'; % RO2 + HO2 = RO
k(:,i) = KRO2HO2.*0.914.*(0.5.*(1.0-fROOH_sec) + 0.5.*(1.0-fROOH_tert));
Gstr{i,1} = 'C9RO2'; Gstr{i,2} = 'HO2'; 
fC9RO2(i)=fC9RO2(i)-1; fHO2(i)=fHO2(i)-1;ffragRO2(i)=ffragRO2(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'C9RO2 + NO = fragRO2 + NO2'; 
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'C9RO2'; Gstr{i,2} = 'NO'; 
fC9RO2(i)=fC9RO2(i)-1; fNO(i)=fNO(i)-1;ffragRO2(i)=ffragRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'C9RO2 + NO = C9NO3'; 
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'C9RO2'; Gstr{i,2} = 'NO'; 
fC9RO2(i)=fC9RO2(i)-1; fNO(i)=fNO(i)-1;fC9NO3(i)=fC9NO3(i)+1;
i=i+1;
Rnames{i} = 'C9RO2 + NO3 = fragRO2 + NO2'; 
k(:,i) = KRO2NO3;
Gstr{i,1} = 'C9RO2'; Gstr{i,2} = 'NO3'; 
fC9RO2(i)=fC9RO2(i)-1; fNO3(i)=fNO3(i)-1;ffragRO2(i)=ffragRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'C9RO2 + RO2 = C9H'; % RO2 + RO2 = ROH
k(:,i) = 2e-12.*(0.5.*bROH_ps + 0.5.*bROH_ta);
Gstr{i,1} = 'C9RO2'; Gstr{i,2} = 'RO2'; 
fC9RO2(i)=fC9RO2(i)-1; fC9H(i)=fC9H(i)+1;
i=i+1;
Rnames{i} = 'C9RO2 + RO2 = C9C'; % RO2 + RO2 = R=O
k(:,i) = 2e-12.*(0.5.*bROH_ps);
Gstr{i,1} = 'C9RO2'; Gstr{i,2} = 'RO2'; 
fC9RO2(i)=fC9RO2(i)-1; fC9C(i)=fC9C(i)+1;
i=i+1;
Rnames{i} = 'C9RO2 + RO2 = fragRO2'; % RO2 + RO2 = RO
k(:,i) = 2e-12.*(0.5.*bRO_ps + 0.5.*bRO_ta);
Gstr{i,1} = 'C9RO2'; Gstr{i,2} = 'RO2'; 
fC9RO2(i)=fC9RO2(i)-1; ffragRO2(i)=ffragRO2(i)+1;

i=i+1;
Rnames{i} = 'C8RO2 + autox = C8RO2';
k(:,i) = kautox_gen; % autoxidation (more oxidized C8RO2)
Gstr{i,1} = 'C8RO2'; 
fC8RO2(i)=fC8RO2(i)-1; fC8RO2(i)=fC8RO2(i)+1;
i=i+1;
Rnames{i} = 'C8RO2 + HO2 = C8OOH'; % RO2 + HO2 = ROOH
k(:,i) = KRO2HO2.*0.914.*(0.5.*fROOH_sec + 0.5.*fROOH_tert);
Gstr{i,1} = 'C8RO2'; Gstr{i,2} = 'HO2'; 
fC8RO2(i)=fC8RO2(i)-1; fHO2(i)=fHO2(i)-1;fC8OOH(i)=fC8OOH(i)+1;
i=i+1;
Rnames{i} = 'C8RO2 + HO2 = fragRO2 + OH'; % RO2 + HO2 = RO
k(:,i) = KRO2HO2.*0.914.*(0.5.*(1.0-fROOH_sec) + 0.5.*(1.0-fROOH_tert));
Gstr{i,1} = 'C8RO2'; Gstr{i,2} = 'HO2'; 
fC8RO2(i)=fC8RO2(i)-1; fHO2(i)=fHO2(i)-1;ffragRO2(i)=ffragRO2(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'C8RO2 + NO = fragRO2 + NO2'; 
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'C8RO2'; Gstr{i,2} = 'NO'; 
fC8RO2(i)=fC8RO2(i)-1; fNO(i)=fNO(i)-1;ffragRO2(i)=ffragRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'C8RO2 + NO = C8NO3'; 
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'C8RO2'; Gstr{i,2} = 'NO'; 
fC8RO2(i)=fC8RO2(i)-1; fNO(i)=fNO(i)-1;fC8NO3(i)=fC8NO3(i)+1;
i=i+1;
Rnames{i} = 'C8RO2 + NO3 = fragRO2 + NO2'; 
k(:,i) = KRO2NO3;
Gstr{i,1} = 'C8RO2'; Gstr{i,2} = 'NO3'; 
fC8RO2(i)=fC8RO2(i)-1; fNO3(i)=fNO3(i)-1;ffragRO2(i)=ffragRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'C8RO2 + RO2 = C8H'; % RO2 + RO2 = ROH
k(:,i) = 2e-12.*(0.5.*bROH_ps + 0.5.*bROH_ta);
Gstr{i,1} = 'C8RO2'; Gstr{i,2} = 'RO2'; 
fC8RO2(i)=fC8RO2(i)-1; fC8H(i)=fC8H(i)+1;
i=i+1;
Rnames{i} = 'C8RO2 + RO2 = C8C'; % RO2 + RO2 = R=O
k(:,i) = 2e-12.*(0.5.*bROH_ps);
Gstr{i,1} = 'C8RO2'; Gstr{i,2} = 'RO2'; 
fC8RO2(i)=fC8RO2(i)-1; fC8C(i)=fC8C(i)+1;
i=i+1;
Rnames{i} = 'C8RO2 + RO2 = fragRO2'; % RO2 + RO2 = RO
k(:,i) = 2e-12.*(0.5.*bRO_ps + 0.5.*bRO_ta);
Gstr{i,1} = 'C8RO2'; Gstr{i,2} = 'RO2'; 
fC8RO2(i)=fC8RO2(i)-1; ffragRO2(i)=ffragRO2(i)+1;

i=i+1;
Rnames{i} = 'C7RO2 + autox = C7RO2';
k(:,i) = kautox_gen; % autoxidation (more oxidized C7RO2)
Gstr{i,1} = 'C7RO2'; 
fC7RO2(i)=fC7RO2(i)-1; fC7RO2(i)=fC7RO2(i)+1;
i=i+1;
Rnames{i} = 'C7RO2 + HO2 = C7OOH'; % RO2 + HO2 = ROOH
k(:,i) = KRO2HO2.*0.914.*fROOH_sec;
Gstr{i,1} = 'C7RO2'; Gstr{i,2} = 'HO2'; 
fC7RO2(i)=fC7RO2(i)-1; fHO2(i)=fHO2(i)-1;fC7OOH(i)=fC7OOH(i)+1;
i=i+1;
Rnames{i} = 'C7RO2 + HO2 = fragRO2 + OH'; % RO2 + HO2 = RO
k(:,i) = KRO2HO2.*0.914.*(1.0-fROOH_sec);
Gstr{i,1} = 'C7RO2'; Gstr{i,2} = 'HO2'; 
fC7RO2(i)=fC7RO2(i)-1; fHO2(i)=fHO2(i)-1;ffragRO2(i)=ffragRO2(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'C7RO2 + NO = fragRO2 + NO2'; 
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'C7RO2'; Gstr{i,2} = 'NO'; 
fC7RO2(i)=fC7RO2(i)-1; fNO(i)=fNO(i)-1;ffragRO2(i)=ffragRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'C7RO2 + NO = C7NO3'; 
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'C7RO2'; Gstr{i,2} = 'NO'; 
fC7RO2(i)=fC7RO2(i)-1; fNO(i)=fNO(i)-1;fC7NO3(i)=fC7NO3(i)+1;
i=i+1;
Rnames{i} = 'C7RO2 + NO3 = fragRO2 + NO2'; 
k(:,i) = KRO2NO3;
Gstr{i,1} = 'C7RO2'; Gstr{i,2} = 'NO3'; 
fC7RO2(i)=fC7RO2(i)-1; fNO3(i)=fNO3(i)-1;ffragRO2(i)=ffragRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'C7RO2 + RO2 = C7H'; % RO2 + RO2 = ROH
k(:,i) = 2e-12.*bROH_ps;
Gstr{i,1} = 'C7RO2'; Gstr{i,2} = 'RO2'; 
fC7RO2(i)=fC7RO2(i)-1; fC7H(i)=fC7H(i)+1;
i=i+1;
Rnames{i} = 'C7RO2 + RO2 = C7C'; % RO2 + RO2 = R=O
k(:,i) = 2e-12.*bROH_ps;
Gstr{i,1} = 'C7RO2'; Gstr{i,2} = 'RO2'; 
fC7RO2(i)=fC7RO2(i)-1; fC7C(i)=fC7C(i)+1;
i=i+1;
Rnames{i} = 'C7RO2 + RO2 = fragRO2'; % RO2 + RO2 = RO
k(:,i) = 2e-12.*bRO_ps;
Gstr{i,1} = 'C7RO2'; Gstr{i,2} = 'RO2'; 
fC7RO2(i)=fC7RO2(i)-1; ffragRO2(i)=ffragRO2(i)+1;

i=i+1;
Rnames{i} = 'fragRO2 + autox = fragRO2';
k(:,i) = kautox_gen; % autoxidation (more oxidized fragRO2)
Gstr{i,1} = 'fragRO2'; 
ffragRO2(i)=ffragRO2(i)-1; ffragRO2(i)=ffragRO2(i)+1;
i=i+1;
Rnames{i} = 'fragRO2 + HO2 = fragOOH'; % RO2 + HO2 = ROOH
k(:,i) = KRO2HO2.*0.914.*(0.5.*fROOH_sec + 0.5.*fROOH_tert);
Gstr{i,1} = 'fragRO2'; Gstr{i,2} = 'HO2'; 
ffragRO2(i)=ffragRO2(i)-1; fHO2(i)=fHO2(i)-1;ffragOOH(i)=ffragOOH(i)+1;
i=i+1;
Rnames{i} = 'fragRO2 + HO2 = fragRO2 + OH'; % RO2 + HO2 = RO
k(:,i) = KRO2HO2.*0.914.*(0.5.*(1.0-fROOH_sec) + 0.5.*(1.0-fROOH_tert));
Gstr{i,1} = 'fragRO2'; Gstr{i,2} = 'HO2'; 
ffragRO2(i)=ffragRO2(i)-1; fHO2(i)=fHO2(i)-1;ffragRO2(i)=ffragRO2(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'fragRO2 + NO = fragRO2 + NO2'; 
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'fragRO2'; Gstr{i,2} = 'NO'; 
ffragRO2(i)=ffragRO2(i)-1; fNO(i)=fNO(i)-1;ffragRO2(i)=ffragRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'fragRO2 + NO = fragNO3'; 
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'fragRO2'; Gstr{i,2} = 'NO'; 
ffragRO2(i)=ffragRO2(i)-1; fNO(i)=fNO(i)-1;ffragNO3(i)=ffragNO3(i)+1;
i=i+1;
Rnames{i} = 'fragRO2 + NO3 = fragRO2 + NO2'; 
k(:,i) = KRO2NO3;
Gstr{i,1} = 'fragRO2'; Gstr{i,2} = 'NO3'; 
ffragRO2(i)=ffragRO2(i)-1; fNO3(i)=fNO3(i)-1;ffragRO2(i)=ffragRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'fragRO2 + RO2 = fragH'; % RO2 + RO2 = ROH
k(:,i) = 2e-12.*(0.5.*bROH_ps + 0.5.*bROH_ta);
Gstr{i,1} = 'fragRO2'; Gstr{i,2} = 'RO2'; 
ffragRO2(i)=ffragRO2(i)-1; ffragH(i)=ffragH(i)+1;
i=i+1;
Rnames{i} = 'fragRO2 + RO2 = fragC'; % RO2 + RO2 = R=O
k(:,i) = 2e-12.*(0.5.*bROH_ps);
Gstr{i,1} = 'fragRO2'; Gstr{i,2} = 'RO2'; 
ffragRO2(i)=ffragRO2(i)-1; ffragC(i)=ffragC(i)+1;
i=i+1;
Rnames{i} = 'fragRO2 + RO2 = fragRO2'; % RO2 + RO2 = RO
k(:,i) = 2e-12.*(0.5.*bRO_ps + 0.5.*bRO_ta);
Gstr{i,1} = 'fragRO2'; Gstr{i,2} = 'RO2'; 
ffragRO2(i)=ffragRO2(i)-1; ffragRO2(i)=ffragRO2(i)+1;

% fragmentation NRO2s (RO2s containing nitrate functional groups)
i=i+1;
Rnames{i} = 'C9NRO2 + autox = C9NRO2';
k(:,i) = kautox_gen; % autoxidation (more oxidized C9RO2)
Gstr{i,1} = 'C9NRO2'; 
fC9NRO2(i)=fC9NRO2(i)-1; fC9NRO2(i)=fC9NRO2(i)+1;
i=i+1;
Rnames{i} = 'C9NRO2 + HO2 = C9NOOH'; % RO2 + HO2 = ROOH
k(:,i) = KRO2HO2.*0.914.*(0.5.*fROOH_sec + 0.5.*fROOH_tert);
Gstr{i,1} = 'C9NRO2'; Gstr{i,2} = 'HO2'; 
fC9NRO2(i)=fC9NRO2(i)-1; fHO2(i)=fHO2(i)-1;fC9NOOH(i)=fC9NOOH(i)+1;
i=i+1;
Rnames{i} = 'C9NRO2 + HO2 = fragNRO2 + OH'; % RO2 + HO2 = RO
k(:,i) = KRO2HO2.*0.914.*(0.5.*(1.0-fROOH_sec) + 0.5.*(1.0-fROOH_tert));
Gstr{i,1} = 'C9NRO2'; Gstr{i,2} = 'HO2'; 
fC9NRO2(i)=fC9NRO2(i)-1; fHO2(i)=fHO2(i)-1;ffragNRO2(i)=ffragNRO2(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'C9NRO2 + NO = fragNRO2 + NO2'; 
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'C9NRO2'; Gstr{i,2} = 'NO'; 
fC9NRO2(i)=fC9NRO2(i)-1; fNO(i)=fNO(i)-1;ffragNRO2(i)=ffragNRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'C9RNO2 + NO = C9NO3'; 
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'C9NRO2'; Gstr{i,2} = 'NO'; 
fC9NRO2(i)=fC9NRO2(i)-1; fNO(i)=fNO(i)-1;fC9NO3(i)=fC9NO3(i)+1;
i=i+1;
Rnames{i} = 'C9NRO2 + NO3 = fragNRO2 + NO2'; 
k(:,i) = KRO2NO3;
Gstr{i,1} = 'C9NRO2'; Gstr{i,2} = 'NO3'; 
fC9NRO2(i)=fC9NRO2(i)-1; fNO3(i)=fNO3(i)-1;ffragNRO2(i)=ffragNRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'C9NRO2 + RO2 = C9NH'; % RO2 + RO2 = ROH
k(:,i) = 2e-12.*(0.5.*bROH_ps + 0.5.*bROH_ta);
Gstr{i,1} = 'C9NRO2'; Gstr{i,2} = 'RO2'; 
fC9NRO2(i)=fC9NRO2(i)-1; fC9NH(i)=fC9NH(i)+1;
i=i+1;
Rnames{i} = 'C9NRO2 + RO2 = C9NC'; % RO2 + RO2 = R=O
k(:,i) = 2e-12.*(0.5.*bROH_ps);
Gstr{i,1} = 'C9NRO2'; Gstr{i,2} = 'RO2'; 
fC9NRO2(i)=fC9NRO2(i)-1; fC9NC(i)=fC9NC(i)+1;
i=i+1;
Rnames{i} = 'C9NRO2 + RO2 = fragRO2'; % RO2 + RO2 = RO
k(:,i) = 2e-12.*(0.5.*bRO_ps + 0.5.*bRO_ta);
Gstr{i,1} = 'C9NRO2'; Gstr{i,2} = 'RO2'; 
fC9NRO2(i)=fC9NRO2(i)-1; ffragNRO2(i)=ffragNRO2(i)+1;

i=i+1;
Rnames{i} = 'C8NRO2 + autox = C8NRO2';
k(:,i) = kautox_gen; % autoxidation (more oxidized C8RO2)
Gstr{i,1} = 'C8NRO2'; 
fC8NRO2(i)=fC8NRO2(i)-1; fC8NRO2(i)=fC8NRO2(i)+1;
i=i+1;
Rnames{i} = 'C8NRO2 + HO2 = C8NOOH'; % RO2 + HO2 = ROOH
k(:,i) = KRO2HO2.*0.914.*(0.5.*fROOH_sec + 0.5.*fROOH_tert);
Gstr{i,1} = 'C8NRO2'; Gstr{i,2} = 'HO2'; 
fC8NRO2(i)=fC8NRO2(i)-1; fHO2(i)=fHO2(i)-1;fC8NOOH(i)=fC8NOOH(i)+1;
i=i+1;
Rnames{i} = 'C8NRO2 + HO2 = fragNRO2 + OH'; % RO2 + HO2 = RO
k(:,i) = KRO2HO2.*0.914.*(0.5.*(1.0-fROOH_sec) + 0.5.*(1.0-fROOH_tert));
Gstr{i,1} = 'C8NRO2'; Gstr{i,2} = 'HO2'; 
fC8NRO2(i)=fC8NRO2(i)-1; fHO2(i)=fHO2(i)-1;ffragNRO2(i)=ffragNRO2(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'C8NRO2 + NO = fragNRO2 + NO2'; 
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'C8NRO2'; Gstr{i,2} = 'NO'; 
fC8NRO2(i)=fC8NRO2(i)-1; fNO(i)=fNO(i)-1;ffragNRO2(i)=ffragNRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'C8NRO2 + NO = C8NO3'; 
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'C8NRO2'; Gstr{i,2} = 'NO'; 
fC8NRO2(i)=fC8NRO2(i)-1; fNO(i)=fNO(i)-1;fC8NO3(i)=fC8NO3(i)+1;
i=i+1;
Rnames{i} = 'C8NRO2 + NO3 = fragNRO2 + NO2'; 
k(:,i) = KRO2NO3;
Gstr{i,1} = 'C8NRO2'; Gstr{i,2} = 'NO3'; 
fC8NRO2(i)=fC8NRO2(i)-1; fNO3(i)=fNO3(i)-1;ffragNRO2(i)=ffragNRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'C8NRO2 + RO2 = C8NH'; % RO2 + RO2 = ROH
k(:,i) = 2e-12.*(0.5.*bROH_ps + 0.5.*bROH_ta);
Gstr{i,1} = 'C8NRO2'; Gstr{i,2} = 'RO2'; 
fC8NRO2(i)=fC8NRO2(i)-1; fC8NH(i)=fC8NH(i)+1;
i=i+1;
Rnames{i} = 'C8NRO2 + RO2 = C8NC'; % RO2 + RO2 = R=O
k(:,i) = 2e-12.*(0.5.*bROH_ps);
Gstr{i,1} = 'C8NRO2'; Gstr{i,2} = 'RO2'; 
fC8NRO2(i)=fC8NRO2(i)-1; fC8NC(i)=fC8NC(i)+1;
i=i+1;
Rnames{i} = 'C8NRO2 + RO2 = fragNRO2'; % RO2 + RO2 = RO
k(:,i) = 2e-12.*(0.5.*bRO_ps + 0.5.*bRO_ta);
Gstr{i,1} = 'C8NRO2'; Gstr{i,2} = 'RO2'; 
fC8NRO2(i)=fC8NRO2(i)-1; ffragNRO2(i)=ffragNRO2(i)+1;

i=i+1;
Rnames{i} = 'C7NRO2 + autox = C7NRO2';
k(:,i) = kautox_gen; % autoxidation (more oxidized C7RO2)
Gstr{i,1} = 'C7NRO2'; 
fC7NRO2(i)=fC7NRO2(i)-1; fC7NRO2(i)=fC7NRO2(i)+1;
i=i+1;
Rnames{i} = 'C7NRO2 + HO2 = C7NOOH'; % RO2 + HO2 = ROOH
k(:,i) = KRO2HO2.*0.914.*(0.5.*fROOH_sec + 0.5.*fROOH_tert);
Gstr{i,1} = 'C7NRO2'; Gstr{i,2} = 'HO2'; 
fC7NRO2(i)=fC7NRO2(i)-1; fHO2(i)=fHO2(i)-1;fC7NOOH(i)=fC7NOOH(i)+1;
i=i+1;
Rnames{i} = 'C7NRO2 + HO2 = fragNRO2 + OH'; % RO2 + HO2 = RO
k(:,i) = KRO2HO2.*0.914.*(0.5.*(1.0-fROOH_sec) + 0.5.*(1.0-fROOH_tert));
Gstr{i,1} = 'C7NRO2'; Gstr{i,2} = 'HO2'; 
fC7NRO2(i)=fC7NRO2(i)-1; fHO2(i)=fHO2(i)-1;ffragNRO2(i)=ffragNRO2(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'C7NRO2 + NO = fragNRO2 + NO2'; 
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'C7NRO2'; Gstr{i,2} = 'NO'; 
fC7NRO2(i)=fC7NRO2(i)-1; fNO(i)=fNO(i)-1;ffragNRO2(i)=ffragNRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'C7NRO2 + NO = C7NO3'; 
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'C7NRO2'; Gstr{i,2} = 'NO'; 
fC7NRO2(i)=fC7NRO2(i)-1; fNO(i)=fNO(i)-1;fC7NO3(i)=fC7NO3(i)+1;
i=i+1;
Rnames{i} = 'C7NRO2 + NO3 = fragNRO2 + NO2'; 
k(:,i) = KRO2NO3;
Gstr{i,1} = 'C7NRO2'; Gstr{i,2} = 'NO3'; 
fC7NRO2(i)=fC7NRO2(i)-1; fNO3(i)=fNO3(i)-1;ffragNRO2(i)=ffragNRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'C7NRO2 + RO2 = C7NH'; % RO2 + RO2 = ROH
k(:,i) = 2e-12.*(0.5.*bROH_ps + 0.5.*bROH_ta);
Gstr{i,1} = 'C7NRO2'; Gstr{i,2} = 'RO2'; 
fC7NRO2(i)=fC7NRO2(i)-1; fC7NH(i)=fC7NH(i)+1;
i=i+1;
Rnames{i} = 'C7NRO2 + RO2 = C7NC'; % RO2 + RO2 = R=O
k(:,i) = 2e-12.*(0.5.*bROH_ps);
Gstr{i,1} = 'C7NRO2'; Gstr{i,2} = 'RO2'; 
fC7NRO2(i)=fC7NRO2(i)-1; fC7NC(i)=fC7NC(i)+1;
i=i+1;
Rnames{i} = 'C7NRO2 + RO2 = fragNRO2'; % RO2 + RO2 = RO
k(:,i) = 2e-12.*(0.5.*bRO_ps + 0.5.*bRO_ta);
Gstr{i,1} = 'C7NRO2'; Gstr{i,2} = 'RO2'; 
fC7NRO2(i)=fC7NRO2(i)-1; ffragNRO2(i)=ffragNRO2(i)+1;

i=i+1;
Rnames{i} = 'fragNRO2 + autox = fragNRO2';
k(:,i) = kautox_gen; % autoxidation (more oxidized fragRO2)
Gstr{i,1} = 'fragNRO2'; 
ffragNRO2(i)=ffragNRO2(i)-1; ffragNRO2(i)=ffragNRO2(i)+1;
i=i+1;
Rnames{i} = 'fragNRO2 + HO2 = fragNOOH'; % RO2 + HO2 = ROOH
k(:,i) = KRO2HO2.*0.914.*(0.5.*fROOH_sec + 0.5.*fROOH_tert);
Gstr{i,1} = 'fragNRO2'; Gstr{i,2} = 'HO2'; 
ffragNRO2(i)=ffragNRO2(i)-1; fHO2(i)=fHO2(i)-1;ffragNOOH(i)=ffragNOOH(i)+1;
i=i+1;
Rnames{i} = 'fragNRO2 + HO2 = fragNRO2 + OH'; % RO2 + HO2 = RO
k(:,i) = KRO2HO2.*0.914.*(0.5.*(1.0-fROOH_sec) + 0.5.*(1.0-fROOH_tert));
Gstr{i,1} = 'fragNRO2'; Gstr{i,2} = 'HO2'; 
ffragNRO2(i)=ffragNRO2(i)-1; fHO2(i)=fHO2(i)-1;ffragNRO2(i)=ffragNRO2(i)+1;fOH(i)=fOH(i)+1;
i=i+1;
Rnames{i} = 'fragNRO2 + NO = fragNRO2 + NO2'; 
k(:,i) = KRO2NO.*0.875;
Gstr{i,1} = 'fragNRO2'; Gstr{i,2} = 'NO'; 
ffragNRO2(i)=ffragNRO2(i)-1; fNO(i)=fNO(i)-1;ffragNRO2(i)=ffragNRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'fragNRO2 + NO = fragNO3'; 
k(:,i) = KRO2NO*0.125;
Gstr{i,1} = 'fragNRO2'; Gstr{i,2} = 'NO'; 
ffragNRO2(i)=ffragNRO2(i)-1; fNO(i)=fNO(i)-1;ffragNO3(i)=ffragNO3(i)+1;
i=i+1;
Rnames{i} = 'fragNRO2 + NO3 = fragNRO2 + NO2'; 
k(:,i) = KRO2NO3;
Gstr{i,1} = 'fragNRO2'; Gstr{i,2} = 'NO3'; 
ffragNRO2(i)=ffragNRO2(i)-1; fNO3(i)=fNO3(i)-1;ffragNRO2(i)=ffragNRO2(i)+1;fNO2(i)=fNO2(i)+1;
i=i+1;
Rnames{i} = 'fragNRO2 + RO2 = fragNH'; % RO2 + RO2 = ROH
k(:,i) = 2e-12.*(0.5.*bROH_ps + 0.5.*bROH_ta);
Gstr{i,1} = 'fragNRO2'; Gstr{i,2} = 'RO2'; 
ffragNRO2(i)=ffragNRO2(i)-1; ffragNH(i)=ffragNH(i)+1;
i=i+1;
Rnames{i} = 'fragNRO2 + RO2 = fragNC'; % RO2 + RO2 = R=O
k(:,i) = 2e-12.*(0.5.*bROH_ps);
Gstr{i,1} = 'fragNRO2'; Gstr{i,2} = 'RO2'; 
ffragNRO2(i)=ffragNRO2(i)-1; ffragNC(i)=ffragNC(i)+1;
i=i+1;
Rnames{i} = 'fragNRO2 + RO2 = fragNRO2'; % RO2 + RO2 = RO
k(:,i) = 2e-12.*(0.5.*bRO_ps + 0.5.*bRO_ta);
Gstr{i,1} = 'fragNRO2'; Gstr{i,2} = 'RO2'; 
ffragNRO2(i)=ffragNRO2(i)-1; ffragNRO2(i)=ffragNRO2(i)+1;

%% dimer formation from RO2 + RO2 = ROOR

%% Dimers fromation between MCM-RO2s
% C107O2 rxns
i=i+1;
Rnames{i} = 'C107O2 + C107O2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C107O2'; 
fC107O2(i)=fC107O2(i)-1; fC107O2(i)=fC107O2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C109O2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C109O2'; 
fC107O2(i)=fC107O2(i)-1; fC109O2(i)=fC109O2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + PINALO2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'PINALO2'; 
fC107O2(i)=fC107O2(i)-1; fPINALO2(i)=fPINALO2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C96CO3 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C96CO3'; 
fC107O2(i)=fC107O2(i)-1; fC96CO3(i)=fC96CO3(i)-1; fC20H30O6(i)=fC20H30O6(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C106O2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C106O2'; 
fC107O2(i)=fC107O2(i)-1; fC106O2(i)=fC106O2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C920CO3 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C920CO3'; 
fC107O2(i)=fC107O2(i)-1; fC920CO3(i)=fC920CO3(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C108O2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C108O2'; 
fC107O2(i)=fC107O2(i)-1; fC108O2(i)=fC108O2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + APINAO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'APINAO2'; 
fC107O2(i)=fC107O2(i)-1; fAPINAO2(i)=fAPINAO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + APINBO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'APINBO2'; 
fC107O2(i)=fC107O2(i)-1; fAPINBO2(i)=fAPINBO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + APINCO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'APINCO2'; 
fC107O2(i)=fC107O2(i)-1; fAPINCO2(i)=fAPINCO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C89CO3 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C89CO3'; 
fC107O2(i)=fC107O2(i)-1; fC89CO3(i)=fC89CO3(i)-1; fC19H28O6(i)=fC19H28O6(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C85CO3 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C85CO3'; 
fC107O2(i)=fC107O2(i)-1; fC85CO3(i)=fC85CO3(i)-1; fC19H28O6(i)=fC19H28O6(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C811CO3 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C811CO3'; 
fC107O2(i)=fC107O2(i)-1; fC811CO3(i)=fC811CO3(i)-1; fC19H28O7(i)=fC19H28O7(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C96O2 = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C96O2'; 
fC107O2(i)=fC107O2(i)-1; fC96O2(i)=fC96O2(i)-1; fC19H30O5(i)=fC19H30O5(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C920O2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C920O2'; 
fC107O2(i)=fC107O2(i)-1; fC920O2(i)=fC920O2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C97O2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C97O2'; 
fC107O2(i)=fC107O2(i)-1; fC97O2(i)=fC97O2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C921O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C921O2'; 
fC107O2(i)=fC107O2(i)-1; fC921O2(i)=fC921O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C98O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C98O2'; 
fC107O2(i)=fC107O2(i)-1; fC98O2(i)=fC98O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C922O2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C922O2'; 
fC107O2(i)=fC107O2(i)-1; fC922O2(i)=fC922O2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C721CO3 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C721CO3'; 
fC107O2(i)=fC107O2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC18H26O7(i)=fC18H26O7(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C85O2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C85O2'; 
fC107O2(i)=fC107O2(i)-1; fC85O2(i)=fC85O2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C89O2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C89O2'; 
fC107O2(i)=fC107O2(i)-1; fC89O2(i)=fC89O2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C86O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C86O2'; 
fC107O2(i)=fC107O2(i)-1; fC86O2(i)=fC86O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C811O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C811O2'; 
fC107O2(i)=fC107O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C810O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C810O2'; 
fC107O2(i)=fC107O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C812O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C812O2'; 
fC107O2(i)=fC107O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + C813O2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'C813O2'; 
fC107O2(i)=fC107O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC107O2(i)=fC107O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C107O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'CHEXO2'; 
fC107O2(i)=fC107O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C109O2 rxns
i=i+1;
Rnames{i} = 'C109O2 + C109O2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C109O2'; 
fC109O2(i)=fC109O2(i)-1; fC109O2(i)=fC109O2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + PINALO2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'PINALO2'; 
fC109O2(i)=fC109O2(i)-1; fPINALO2(i)=fPINALO2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C96CO3 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C96CO3'; 
fC109O2(i)=fC109O2(i)-1; fC96CO3(i)=fC96CO3(i)-1; fC20H30O6(i)=fC20H30O6(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C106O2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C106O2'; 
fC109O2(i)=fC109O2(i)-1; fC106O2(i)=fC106O2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C920CO3 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C920CO3'; 
fC109O2(i)=fC109O2(i)-1; fC920CO3(i)=fC920CO3(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C108O2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C108O2'; 
fC109O2(i)=fC109O2(i)-1; fC108O2(i)=fC108O2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + APINAO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'APINAO2'; 
fC109O2(i)=fC109O2(i)-1; fAPINAO2(i)=fAPINAO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + APINBO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'APINBO2'; 
fC109O2(i)=fC109O2(i)-1; fAPINBO2(i)=fAPINBO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + APINCO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'APINCO2'; 
fC109O2(i)=fC109O2(i)-1; fAPINCO2(i)=fAPINCO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C89CO3 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C89CO3'; 
fC109O2(i)=fC109O2(i)-1; fC89CO3(i)=fC89CO3(i)-1; fC19H28O6(i)=fC19H28O6(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C85CO3 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C85CO3'; 
fC109O2(i)=fC109O2(i)-1; fC85CO3(i)=fC85CO3(i)-1; fC19H28O6(i)=fC19H28O6(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C811CO3 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C811CO3'; 
fC109O2(i)=fC109O2(i)-1; fC811CO3(i)=fC811CO3(i)-1; fC19H28O7(i)=fC19H28O7(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C96O2 = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C96O2'; 
fC109O2(i)=fC109O2(i)-1; fC96O2(i)=fC96O2(i)-1; fC19H30O5(i)=fC19H30O5(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C920O2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C920O2'; 
fC109O2(i)=fC109O2(i)-1; fC920O2(i)=fC920O2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C97O2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C97O2'; 
fC109O2(i)=fC109O2(i)-1; fC97O2(i)=fC97O2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C921O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C921O2'; 
fC109O2(i)=fC109O2(i)-1; fC921O2(i)=fC921O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C98O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C98O2'; 
fC109O2(i)=fC109O2(i)-1; fC98O2(i)=fC98O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C922O2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C922O2'; 
fC109O2(i)=fC109O2(i)-1; fC922O2(i)=fC922O2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C721CO3 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C721CO3'; 
fC109O2(i)=fC109O2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC18H26O7(i)=fC18H26O7(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C85O2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C85O2'; 
fC109O2(i)=fC109O2(i)-1; fC85O2(i)=fC85O2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C89O2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C89O2'; 
fC109O2(i)=fC109O2(i)-1; fC89O2(i)=fC89O2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C86O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C86O2'; 
fC109O2(i)=fC109O2(i)-1; fC86O2(i)=fC86O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C811O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C811O2'; 
fC109O2(i)=fC109O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C810O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C810O2'; 
fC109O2(i)=fC109O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C812O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C812O2'; 
fC109O2(i)=fC109O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + C813O2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'C813O2'; 
fC109O2(i)=fC109O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC109O2(i)=fC109O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C109O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'CHEXO2'; 
fC109O2(i)=fC109O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% PINALO2 rxns
i=i+1;
Rnames{i} = 'PINALO2 + PINALO2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'PINALO2'; 
fPINALO2(i)=fPINALO2(i)-1; fPINALO2(i)=fPINALO2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C96CO3 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C96CO3'; 
fPINALO2(i)=fPINALO2(i)-1; fC96CO3(i)=fC96CO3(i)-1; fC20H30O6(i)=fC20H30O6(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C106O2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C106O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC106O2(i)=fC106O2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C920CO3 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C920CO3'; 
fPINALO2(i)=fPINALO2(i)-1; fC920CO3(i)=fC920CO3(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C108O2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C108O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC108O2(i)=fC108O2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + APINAO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'APINAO2'; 
fPINALO2(i)=fPINALO2(i)-1; fAPINAO2(i)=fAPINAO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + APINBO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'APINBO2'; 
fPINALO2(i)=fPINALO2(i)-1; fAPINBO2(i)=fAPINBO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + APINCO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'APINCO2'; 
fPINALO2(i)=fPINALO2(i)-1; fAPINCO2(i)=fAPINCO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C89CO3 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C89CO3'; 
fPINALO2(i)=fPINALO2(i)-1; fC89CO3(i)=fC89CO3(i)-1; fC19H28O6(i)=fC19H28O6(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C85CO3 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C85CO3'; 
fPINALO2(i)=fPINALO2(i)-1; fC85CO3(i)=fC85CO3(i)-1; fC19H28O6(i)=fC19H28O6(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C811CO3 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C811CO3'; 
fPINALO2(i)=fPINALO2(i)-1; fC811CO3(i)=fC811CO3(i)-1; fC19H28O7(i)=fC19H28O7(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C96O2 = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C96O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC96O2(i)=fC96O2(i)-1; fC19H30O5(i)=fC19H30O5(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C920O2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C920O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC920O2(i)=fC920O2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C97O2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C97O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC97O2(i)=fC97O2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C921O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C921O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC921O2(i)=fC921O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C98O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C98O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC98O2(i)=fC98O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C922O2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C922O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC922O2(i)=fC922O2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C721CO3 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C721CO3'; 
fPINALO2(i)=fPINALO2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC18H26O7(i)=fC18H26O7(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C85O2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C85O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC85O2(i)=fC85O2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C89O2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C89O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC89O2(i)=fC89O2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C86O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C86O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC86O2(i)=fC86O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C811O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C811O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC811O2(i)=fC811O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C810O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C810O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC810O2(i)=fC810O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C812O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C812O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC812O2(i)=fC812O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + C813O2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'C813O2'; 
fPINALO2(i)=fPINALO2(i)-1; fC813O2(i)=fC813O2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'BUT2OLO2'; 
fPINALO2(i)=fPINALO2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'PINALO2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'CHEXO2'; 
fPINALO2(i)=fPINALO2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C96CO3 rxns
i=i+1;
Rnames{i} = 'C96CO3 + C96CO3 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C96CO3'; 
fC96CO3(i)=fC96CO3(i)-1; fC96CO3(i)=fC96CO3(i)-1; fC20H30O6(i)=fC20H30O6(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C106O2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C106O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC106O2(i)=fC106O2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C920CO3 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C920CO3'; 
fC96CO3(i)=fC96CO3(i)-1; fC920CO3(i)=fC920CO3(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C108O2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C108O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC108O2(i)=fC108O2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + APINAO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'APINAO2'; 
fC96CO3(i)=fC96CO3(i)-1; fAPINAO2(i)=fAPINAO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + APINBO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'APINBO2'; 
fC96CO3(i)=fC96CO3(i)-1; fAPINBO2(i)=fAPINBO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + APINCO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'APINCO2'; 
fC96CO3(i)=fC96CO3(i)-1; fAPINCO2(i)=fAPINCO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C89CO3 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C89CO3'; 
fC96CO3(i)=fC96CO3(i)-1; fC89CO3(i)=fC89CO3(i)-1; fC19H28O6(i)=fC19H28O6(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C85CO3 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C85CO3'; 
fC96CO3(i)=fC96CO3(i)-1; fC85CO3(i)=fC85CO3(i)-1; fC19H28O6(i)=fC19H28O6(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C811CO3 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C811CO3'; 
fC96CO3(i)=fC96CO3(i)-1; fC811CO3(i)=fC811CO3(i)-1; fC19H28O7(i)=fC19H28O7(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C96O2 = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C96O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC96O2(i)=fC96O2(i)-1; fC19H30O5(i)=fC19H30O5(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C920O2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C920O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC920O2(i)=fC920O2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C97O2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C97O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC97O2(i)=fC97O2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C921O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C921O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC921O2(i)=fC921O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C98O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C98O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC98O2(i)=fC98O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C922O2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C922O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC922O2(i)=fC922O2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C721CO3 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C721CO3'; 
fC96CO3(i)=fC96CO3(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC18H26O7(i)=fC18H26O7(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C85O2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C85O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC85O2(i)=fC85O2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C89O2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C89O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC89O2(i)=fC89O2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C86O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C86O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC86O2(i)=fC86O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C811O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C811O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC811O2(i)=fC811O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C810O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C810O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC810O2(i)=fC810O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C812O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C812O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC812O2(i)=fC812O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + C813O2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'C813O2'; 
fC96CO3(i)=fC96CO3(i)-1; fC813O2(i)=fC813O2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'BUT2OLO2'; 
fC96CO3(i)=fC96CO3(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C96CO3 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'CHEXO2'; 
fC96CO3(i)=fC96CO3(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;


% C106O2 rxns
i=i+1;
Rnames{i} = 'C106O2 + C106O2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C106O2'; 
fC106O2(i)=fC106O2(i)-1; fC106O2(i)=fC106O2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C920CO3 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C920CO3'; 
fC106O2(i)=fC106O2(i)-1; fC920CO3(i)=fC920CO3(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C108O2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C108O2'; 
fC106O2(i)=fC106O2(i)-1; fC108O2(i)=fC108O2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + APINAO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'APINAO2'; 
fC106O2(i)=fC106O2(i)-1; fAPINAO2(i)=fAPINAO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + APINBO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'APINBO2'; 
fC106O2(i)=fC106O2(i)-1; fAPINBO2(i)=fAPINBO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + APINCO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'APINCO2'; 
fC106O2(i)=fC106O2(i)-1; fAPINCO2(i)=fAPINCO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C89CO3 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C89CO3'; 
fC106O2(i)=fC106O2(i)-1; fC89CO3(i)=fC89CO3(i)-1; fC19H28O7(i)=fC19H28O7(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C85CO3 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C85CO3'; 
fC106O2(i)=fC106O2(i)-1; fC85CO3(i)=fC85CO3(i)-1; fC19H28O7(i)=fC19H28O7(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C811CO3 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C811CO3'; 
fC106O2(i)=fC106O2(i)-1; fC811CO3(i)=fC811CO3(i)-1; fC19H28O8(i)=fC19H28O8(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C96O2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C96O2'; 
fC106O2(i)=fC106O2(i)-1; fC96O2(i)=fC96O2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C920O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C920O2'; 
fC106O2(i)=fC106O2(i)-1; fC920O2(i)=fC920O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C97O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C97O2'; 
fC106O2(i)=fC106O2(i)-1; fC97O2(i)=fC97O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C921O2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C921O2'; 
fC106O2(i)=fC106O2(i)-1; fC921O2(i)=fC921O2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C98O2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C98O2'; 
fC106O2(i)=fC106O2(i)-1; fC98O2(i)=fC98O2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C922O2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C922O2'; 
fC106O2(i)=fC106O2(i)-1; fC922O2(i)=fC922O2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C721CO3 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C721CO3'; 
fC106O2(i)=fC106O2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC18H26O8(i)=fC18H26O8(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C85O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C85O2'; 
fC106O2(i)=fC106O2(i)-1; fC85O2(i)=fC85O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C89O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C89O2'; 
fC106O2(i)=fC106O2(i)-1; fC89O2(i)=fC89O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C86O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C86O2'; 
fC106O2(i)=fC106O2(i)-1; fC86O2(i)=fC86O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C811O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C811O2'; 
fC106O2(i)=fC106O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C810O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C810O2'; 
fC106O2(i)=fC106O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C812O2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C812O2'; 
fC106O2(i)=fC106O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + C813O2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'C813O2'; 
fC106O2(i)=fC106O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC106O2(i)=fC106O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C106O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'CHEXO2'; 
fC106O2(i)=fC106O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C920CO3 rxns
i=i+1;
Rnames{i} = 'C920CO3 + C920CO3 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C920CO3'; 
fC920CO3(i)=fC920CO3(i)-1; fC920CO3(i)=fC920CO3(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C108O2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C108O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC108O2(i)=fC108O2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + APINAO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'APINAO2'; 
fC920CO3(i)=fC920CO3(i)-1; fAPINAO2(i)=fAPINAO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + APINBO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'APINBO2'; 
fC920CO3(i)=fC920CO3(i)-1; fAPINBO2(i)=fAPINBO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + APINCO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'APINCO2'; 
fC920CO3(i)=fC920CO3(i)-1; fAPINCO2(i)=fAPINCO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C89CO3 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C89CO3'; 
fC920CO3(i)=fC920CO3(i)-1; fC89CO3(i)=fC89CO3(i)-1; fC19H28O7(i)=fC19H28O7(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C85CO3 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C85CO3'; 
fC920CO3(i)=fC920CO3(i)-1; fC85CO3(i)=fC85CO3(i)-1; fC19H28O7(i)=fC19H28O7(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C811CO3 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C811CO3'; 
fC920CO3(i)=fC920CO3(i)-1; fC811CO3(i)=fC811CO3(i)-1; fC19H28O8(i)=fC19H28O8(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C96O2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C96O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC96O2(i)=fC96O2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C920O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C920O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC920O2(i)=fC920O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C97O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C97O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC97O2(i)=fC97O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C921O2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C921O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC921O2(i)=fC921O2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C98O2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C98O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC98O2(i)=fC98O2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C922O2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C922O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC922O2(i)=fC922O2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C721CO3 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C721CO3'; 
fC920CO3(i)=fC920CO3(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC18H26O8(i)=fC18H26O8(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C85O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C85O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC85O2(i)=fC85O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C89O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C89O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC89O2(i)=fC89O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C86O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C86O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC86O2(i)=fC86O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C811O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C811O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC811O2(i)=fC811O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C810O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C810O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC810O2(i)=fC810O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C812O2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C812O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC812O2(i)=fC812O2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + C813O2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'C813O2'; 
fC920CO3(i)=fC920CO3(i)-1; fC813O2(i)=fC813O2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'BUT2OLO2'; 
fC920CO3(i)=fC920CO3(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C920CO3 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'CHEXO2'; 
fC920CO3(i)=fC920CO3(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C108O2 rxns
i=i+1;
Rnames{i} = 'C108O2 + C108O2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C108O2'; 
fC108O2(i)=fC108O2(i)-1; fC108O2(i)=fC108O2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + APINAO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'APINAO2'; 
fC108O2(i)=fC108O2(i)-1; fAPINAO2(i)=fAPINAO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + APINBO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'APINBO2'; 
fC108O2(i)=fC108O2(i)-1; fAPINBO2(i)=fAPINBO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + APINCO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'APINCO2'; 
fC108O2(i)=fC108O2(i)-1; fAPINCO2(i)=fAPINCO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C89CO3 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C89CO3'; 
fC108O2(i)=fC108O2(i)-1; fC89CO3(i)=fC89CO3(i)-1; fC19H28O7(i)=fC19H28O7(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C85CO3 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C85CO3'; 
fC108O2(i)=fC108O2(i)-1; fC85CO3(i)=fC85CO3(i)-1; fC19H28O7(i)=fC19H28O7(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C811CO3 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C811CO3'; 
fC108O2(i)=fC108O2(i)-1; fC811CO3(i)=fC811CO3(i)-1; fC19H28O8(i)=fC19H28O8(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C96O2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C96O2'; 
fC108O2(i)=fC108O2(i)-1; fC96O2(i)=fC96O2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C920O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C920O2'; 
fC108O2(i)=fC108O2(i)-1; fC920O2(i)=fC920O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C97O2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C97O2'; 
fC108O2(i)=fC108O2(i)-1; fC97O2(i)=fC97O2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C921O2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C921O2'; 
fC108O2(i)=fC108O2(i)-1; fC921O2(i)=fC921O2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C98O2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C98O2'; 
fC108O2(i)=fC108O2(i)-1; fC98O2(i)=fC98O2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C922O2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C922O2'; 
fC108O2(i)=fC108O2(i)-1; fC922O2(i)=fC922O2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C721CO3 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C721CO3'; 
fC108O2(i)=fC108O2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC18H26O8(i)=fC18H26O8(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C85O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C85O2'; 
fC108O2(i)=fC108O2(i)-1; fC85O2(i)=fC85O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C89O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C89O2'; 
fC108O2(i)=fC108O2(i)-1; fC89O2(i)=fC89O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C86O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C86O2'; 
fC108O2(i)=fC108O2(i)-1; fC86O2(i)=fC86O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C811O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C811O2'; 
fC108O2(i)=fC108O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C810O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C810O2'; 
fC108O2(i)=fC108O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C812O2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C812O2'; 
fC108O2(i)=fC108O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + C813O2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'C813O2'; 
fC108O2(i)=fC108O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC108O2(i)=fC108O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C108O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'CHEXO2'; 
fC108O2(i)=fC108O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% APINAO2 rxns
i=i+1;
Rnames{i} = 'APINAO2 + APINAO2 = C20H34O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'APINAO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fAPINAO2(i)=fAPINAO2(i)-1; fC20H34O4(i)=fC20H34O4(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + APINBO2 = C20H34O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'APINBO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fAPINBO2(i)=fAPINBO2(i)-1; fC20H34O4(i)=fC20H34O4(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + APINCO2 = C20H34O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'APINCO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fAPINCO2(i)=fAPINCO2(i)-1; fC20H34O4(i)=fC20H34O4(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C89CO3 = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C89CO3'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC89CO3(i)=fC89CO3(i)-1; fC19H30O5(i)=fC19H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C85CO3 = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C85CO3'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC85CO3(i)=fC85CO3(i)-1; fC19H30O5(i)=fC19H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C811CO3 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C811CO3'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC811CO3(i)=fC811CO3(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C96O2 = C19H32O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C96O2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC96O2(i)=fC96O2(i)-1; fC19H32O4(i)=fC19H32O4(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C920O2 = C19H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C920O2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC920O2(i)=fC920O2(i)-1; fC19H32O5(i)=fC19H32O5(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C97O2 = C19H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C97O2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC97O2(i)=fC97O2(i)-1; fC19H32O5(i)=fC19H32O5(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C921O2 = C19H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C921O2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC921O2(i)=fC921O2(i)-1; fC19H32O6(i)=fC19H32O6(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C98O2 = C19H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C98O2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC98O2(i)=fC98O2(i)-1; fC19H32O6(i)=fC19H32O6(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C922O2 = C19H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C922O2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC922O2(i)=fC922O2(i)-1; fC19H32O7(i)=fC19H32O7(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C721CO3 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C721CO3'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C85O2 = C18H30O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C85O2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC85O2(i)=fC85O2(i)-1; fC18H30O4(i)=fC18H30O4(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C89O2 = C18H30O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C89O2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC89O2(i)=fC89O2(i)-1; fC18H30O4(i)=fC18H30O4(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C86O2 = C18H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C86O2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC86O2(i)=fC86O2(i)-1; fC18H30O5(i)=fC18H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C811O2 = C18H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C811O2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC811O2(i)=fC811O2(i)-1; fC18H30O5(i)=fC18H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C810O2 = C18H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C810O2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC810O2(i)=fC810O2(i)-1; fC18H30O5(i)=fC18H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C812O2 = C18H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C812O2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC812O2(i)=fC812O2(i)-1; fC18H30O6(i)=fC18H30O6(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + C813O2 = C18H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'C813O2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fC813O2(i)=fC813O2(i)-1; fC18H30O7(i)=fC18H30O7(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'BUT2OLO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'APINAO2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'CHEXO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% APINBO2 rxns
i=i+1;
Rnames{i} = 'APINBO2 + APINBO2 = C20H34O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'APINBO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fAPINBO2(i)=fAPINBO2(i)-1; fC20H34O4(i)=fC20H34O4(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + APINCO2 = C20H34O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'APINCO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fAPINCO2(i)=fAPINCO2(i)-1; fC20H34O4(i)=fC20H34O4(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C89CO3 = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C89CO3'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC89CO3(i)=fC89CO3(i)-1; fC19H30O5(i)=fC19H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C85CO3 = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C85CO3'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC85CO3(i)=fC85CO3(i)-1; fC19H30O5(i)=fC19H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C811CO3 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C811CO3'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC811CO3(i)=fC811CO3(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C96O2 = C19H32O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C96O2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC96O2(i)=fC96O2(i)-1; fC19H32O4(i)=fC19H32O4(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C920O2 = C19H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C920O2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC920O2(i)=fC920O2(i)-1; fC19H32O5(i)=fC19H32O5(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C97O2 = C19H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C97O2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC97O2(i)=fC97O2(i)-1; fC19H32O5(i)=fC19H32O5(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C921O2 = C19H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C921O2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC921O2(i)=fC921O2(i)-1; fC19H32O6(i)=fC19H32O6(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C98O2 = C19H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C98O2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC98O2(i)=fC98O2(i)-1; fC19H32O6(i)=fC19H32O6(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C922O2 = C19H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C922O2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC922O2(i)=fC922O2(i)-1; fC19H32O7(i)=fC19H32O7(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C721CO3 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C721CO3'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C85O2 = C18H30O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C85O2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC85O2(i)=fC85O2(i)-1; fC18H30O4(i)=fC18H30O4(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C89O2 = C18H30O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C89O2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC89O2(i)=fC89O2(i)-1; fC18H30O4(i)=fC18H30O4(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C86O2 = C18H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C86O2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC86O2(i)=fC86O2(i)-1; fC18H30O5(i)=fC18H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C811O2 = C18H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C811O2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC811O2(i)=fC811O2(i)-1; fC18H30O5(i)=fC18H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C810O2 = C18H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C810O2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC810O2(i)=fC810O2(i)-1; fC18H30O5(i)=fC18H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C812O2 = C18H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C812O2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC812O2(i)=fC812O2(i)-1; fC18H30O6(i)=fC18H30O6(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + C813O2 = C18H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'C813O2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fC813O2(i)=fC813O2(i)-1; fC18H30O7(i)=fC18H30O7(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'BUT2OLO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'APINBO2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'CHEXO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% APINCO2 rxns
i=i+1;
Rnames{i} = 'APINCO2 + APINCO2 = C20H34O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'APINCO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fAPINCO2(i)=fAPINCO2(i)-1; fC20H34O4(i)=fC20H34O4(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C89CO3 = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C89CO3'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC89CO3(i)=fC89CO3(i)-1; fC19H30O5(i)=fC19H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C85CO3 = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C85CO3'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC85CO3(i)=fC85CO3(i)-1; fC19H30O5(i)=fC19H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C811CO3 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C811CO3'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC811CO3(i)=fC811CO3(i)-1; fC19H30O6(i)=fC19H30O6(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C96O2 = C19H32O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C96O2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC96O2(i)=fC96O2(i)-1; fC19H32O4(i)=fC19H32O4(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C920O2 = C19H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C920O2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC920O2(i)=fC920O2(i)-1; fC19H32O5(i)=fC19H32O5(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C97O2 = C19H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C97O2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC97O2(i)=fC97O2(i)-1; fC19H32O5(i)=fC19H32O5(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C921O2 = C19H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C921O2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC921O2(i)=fC921O2(i)-1; fC19H32O6(i)=fC19H32O6(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C98O2 = C19H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C98O2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC98O2(i)=fC98O2(i)-1; fC19H32O6(i)=fC19H32O6(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C922O2 = C19H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C922O2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC922O2(i)=fC922O2(i)-1; fC19H32O7(i)=fC19H32O7(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C721CO3 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C721CO3'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C85O2 = C18H30O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C85O2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC85O2(i)=fC85O2(i)-1; fC18H30O4(i)=fC18H30O4(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C89O2 = C18H30O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C89O2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC89O2(i)=fC89O2(i)-1; fC18H30O4(i)=fC18H30O4(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C86O2 = C18H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C86O2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC86O2(i)=fC86O2(i)-1; fC18H30O5(i)=fC18H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C811O2 = C18H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C811O2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC811O2(i)=fC811O2(i)-1; fC18H30O5(i)=fC18H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C810O2 = C18H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C810O2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC810O2(i)=fC810O2(i)-1; fC18H30O5(i)=fC18H30O5(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C812O2 = C18H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C812O2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC812O2(i)=fC812O2(i)-1; fC18H30O6(i)=fC18H30O6(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + C813O2 = C18H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'C813O2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fC813O2(i)=fC813O2(i)-1; fC18H30O7(i)=fC18H30O7(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'BUT2OLO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'APINCO2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'CHEXO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C89CO3 rxns
i=i+1;
Rnames{i} = 'C89CO3 + C89CO3 = C18H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C89CO3'; 
fC89CO3(i)=fC89CO3(i)-1; fC89CO3(i)=fC89CO3(i)-1; fC18H26O6(i)=fC18H26O6(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C85CO3 = C18H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C85CO3'; 
fC89CO3(i)=fC89CO3(i)-1; fC85CO3(i)=fC85CO3(i)-1; fC18H26O6(i)=fC18H26O6(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C811CO3 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C811CO3'; 
fC89CO3(i)=fC89CO3(i)-1; fC811CO3(i)=fC811CO3(i)-1; fC18H26O7(i)=fC18H26O7(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C96O2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C96O2'; 
fC89CO3(i)=fC89CO3(i)-1; fC96O2(i)=fC96O2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C920O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C920O2'; 
fC89CO3(i)=fC89CO3(i)-1; fC920O2(i)=fC920O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C97O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C97O2'; 
fC89CO3(i)=fC89CO3(i)-1; fC97O2(i)=fC97O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C921O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C921O2'; 
fC89CO3(i)=fC89CO3(i)-1; fC921O2(i)=fC921O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C98O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C98O2'; 
fC89CO3(i)=fC89CO3(i)-1; fC98O2(i)=fC98O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C922O2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C922O2'; 
fC89CO3(i)=fC89CO3(i)-1; fC922O2(i)=fC922O2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C721CO3 = C17H24O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C721CO3'; 
fC89CO3(i)=fC89CO3(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC17H24O7(i)=fC17H24O7(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C85O2 = C17H26O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C85O2'; 
fC89CO3(i)=fC89CO3(i)-1; fC85O2(i)=fC85O2(i)-1; fC17H26O5(i)=fC17H26O5(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C89O2 = C17H26O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C89O2'; 
fC89CO3(i)=fC89CO3(i)-1; fC89O2(i)=fC89O2(i)-1; fC17H26O5(i)=fC17H26O5(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C86O2 = C17H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C86O2'; 
fC89CO3(i)=fC89CO3(i)-1; fC86O2(i)=fC86O2(i)-1; fC17H26O6(i)=fC17H26O6(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C811O2 = C17H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C811O2'; 
fC89CO3(i)=fC89CO3(i)-1; fC811O2(i)=fC811O2(i)-1; fC17H26O6(i)=fC17H26O6(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C810O2 = C17H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C810O2'; 
fC89CO3(i)=fC89CO3(i)-1; fC810O2(i)=fC810O2(i)-1; fC17H26O6(i)=fC17H26O6(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C812O2 = C17H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C812O2'; 
fC89CO3(i)=fC89CO3(i)-1; fC812O2(i)=fC812O2(i)-1; fC17H26O7(i)=fC17H26O7(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + C813O2 = C17H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'C813O2'; 
fC89CO3(i)=fC89CO3(i)-1; fC813O2(i)=fC813O2(i)-1; fC17H26O8(i)=fC17H26O8(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'BUT2OLO2'; 
fC89CO3(i)=fC89CO3(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C89CO3 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'CHEXO2'; 
fC89CO3(i)=fC89CO3(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C85CO3 rxns
i=i+1;
Rnames{i} = 'C85CO3 + C85CO3 = C18H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C85CO3'; 
fC85CO3(i)=fC85CO3(i)-1; fC85CO3(i)=fC85CO3(i)-1; fC18H26O6(i)=fC18H26O6(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C811CO3 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C811CO3'; 
fC85CO3(i)=fC85CO3(i)-1; fC811CO3(i)=fC811CO3(i)-1; fC18H26O7(i)=fC18H26O7(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C96O2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C96O2'; 
fC85CO3(i)=fC85CO3(i)-1; fC96O2(i)=fC96O2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C920O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C920O2'; 
fC85CO3(i)=fC85CO3(i)-1; fC920O2(i)=fC920O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C97O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C97O2'; 
fC85CO3(i)=fC85CO3(i)-1; fC97O2(i)=fC97O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C921O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C921O2'; 
fC85CO3(i)=fC85CO3(i)-1; fC921O2(i)=fC921O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C98O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C98O2'; 
fC85CO3(i)=fC85CO3(i)-1; fC98O2(i)=fC98O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C922O2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C922O2'; 
fC85CO3(i)=fC85CO3(i)-1; fC922O2(i)=fC922O2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C721CO3 = C17H24O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C721CO3'; 
fC85CO3(i)=fC85CO3(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC17H24O7(i)=fC17H24O7(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C85O2 = C17H26O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C85O2'; 
fC85CO3(i)=fC85CO3(i)-1; fC85O2(i)=fC85O2(i)-1; fC17H26O5(i)=fC17H26O5(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C89O2 = C17H26O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C89O2'; 
fC85CO3(i)=fC85CO3(i)-1; fC89O2(i)=fC89O2(i)-1; fC17H26O5(i)=fC17H26O5(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C86O2 = C17H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C86O2'; 
fC85CO3(i)=fC85CO3(i)-1; fC86O2(i)=fC86O2(i)-1; fC17H26O6(i)=fC17H26O6(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C811O2 = C17H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C811O2'; 
fC85CO3(i)=fC85CO3(i)-1; fC811O2(i)=fC811O2(i)-1; fC17H26O6(i)=fC17H26O6(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C810O2 = C17H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C810O2'; 
fC85CO3(i)=fC85CO3(i)-1; fC810O2(i)=fC810O2(i)-1; fC17H26O6(i)=fC17H26O6(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C812O2 = C17H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C812O2'; 
fC85CO3(i)=fC85CO3(i)-1; fC812O2(i)=fC812O2(i)-1; fC17H26O7(i)=fC17H26O7(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + C813O2 = C17H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'C813O2'; 
fC85CO3(i)=fC85CO3(i)-1; fC813O2(i)=fC813O2(i)-1; fC17H26O8(i)=fC17H26O8(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'BUT2OLO2'; 
fC85CO3(i)=fC85CO3(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C85CO3 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'CHEXO2'; 
fC85CO3(i)=fC85CO3(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C811CO3
i=i+1;
Rnames{i} = 'C811CO3 + C811CO3 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C811CO3'; 
fC811CO3(i)=fC811CO3(i)-1; fC811CO3(i)=fC811CO3(i)-1; fC18H26O8(i)=fC18H26O8(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C96O2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C96O2'; 
fC811CO3(i)=fC811CO3(i)-1; fC96O2(i)=fC96O2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C920O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C920O2'; 
fC811CO3(i)=fC811CO3(i)-1; fC920O2(i)=fC920O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C97O2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C97O2'; 
fC811CO3(i)=fC811CO3(i)-1; fC97O2(i)=fC97O2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C921O2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C921O2'; 
fC811CO3(i)=fC811CO3(i)-1; fC921O2(i)=fC921O2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C98O2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C98O2'; 
fC811CO3(i)=fC811CO3(i)-1; fC98O2(i)=fC98O2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C922O2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C922O2'; 
fC811CO3(i)=fC811CO3(i)-1; fC922O2(i)=fC922O2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C721CO3 = C17H24O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C721CO3'; 
fC811CO3(i)=fC811CO3(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC17H24O8(i)=fC17H24O8(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C85O2 = C17H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C85O2'; 
fC811CO3(i)=fC811CO3(i)-1; fC85O2(i)=fC85O2(i)-1; fC17H26O6(i)=fC17H26O6(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C89O2 = C17H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C89O2'; 
fC811CO3(i)=fC811CO3(i)-1; fC89O2(i)=fC89O2(i)-1; fC17H26O6(i)=fC17H26O6(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C86O2 = C17H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C86O2'; 
fC811CO3(i)=fC811CO3(i)-1; fC86O2(i)=fC86O2(i)-1; fC17H26O7(i)=fC17H26O7(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C811O2 = C17H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C811O2'; 
fC811CO3(i)=fC811CO3(i)-1; fC811O2(i)=fC811O2(i)-1; fC17H26O7(i)=fC17H26O7(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C810O2 = C17H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C810O2'; 
fC811CO3(i)=fC811CO3(i)-1; fC810O2(i)=fC810O2(i)-1; fC17H26O7(i)=fC17H26O7(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C812O2 = C17H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C812O2'; 
fC811CO3(i)=fC811CO3(i)-1; fC812O2(i)=fC812O2(i)-1; fC17H26O8(i)=fC17H26O8(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + C813O2 = C17H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'C813O2'; 
fC811CO3(i)=fC811CO3(i)-1; fC813O2(i)=fC813O2(i)-1; fC17H26O9(i)=fC17H26O9(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'BUT2OLO2'; 
fC811CO3(i)=fC811CO3(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C811CO3 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'CHEXO2'; 
fC811CO3(i)=fC811CO3(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C96O2 rxns
i=i+1;
Rnames{i} = 'C96O2 + C96O2 = C18H30O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C96O2'; 
fC96O2(i)=fC96O2(i)-1; fC96O2(i)=fC96O2(i)-1; fC18H30O4(i)=fC18H30O4(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + C920O2 = C18H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C920O2'; 
fC96O2(i)=fC96O2(i)-1; fC920O2(i)=fC920O2(i)-1; fC18H30O5(i)=fC18H30O5(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + C97O2 = C18H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C97O2'; 
fC96O2(i)=fC96O2(i)-1; fC97O2(i)=fC97O2(i)-1; fC18H30O5(i)=fC18H30O5(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + C921O2 = C18H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C921O2'; 
fC96O2(i)=fC96O2(i)-1; fC921O2(i)=fC921O2(i)-1; fC18H30O6(i)=fC18H30O6(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + C98O2 = C18H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C98O2'; 
fC96O2(i)=fC96O2(i)-1; fC98O2(i)=fC98O2(i)-1; fC18H30O6(i)=fC18H30O6(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + C922O2 = C18H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C922O2'; 
fC96O2(i)=fC96O2(i)-1; fC922O2(i)=fC922O2(i)-1; fC18H30O7(i)=fC18H30O7(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + C721CO3 = C17H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C721CO3'; 
fC96O2(i)=fC96O2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC17H26O6(i)=fC17H26O6(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + C85O2 = C17H28O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C85O2'; 
fC96O2(i)=fC96O2(i)-1; fC85O2(i)=fC85O2(i)-1; fC17H28O4(i)=fC17H28O4(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + C89O2 = C17H28O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C89O2'; 
fC96O2(i)=fC96O2(i)-1; fC89O2(i)=fC89O2(i)-1; fC17H28O4(i)=fC17H28O4(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + C86O2 = C17H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C86O2'; 
fC96O2(i)=fC96O2(i)-1; fC86O2(i)=fC86O2(i)-1; fC17H28O5(i)=fC17H28O5(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + C811O2 = C17H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C811O2'; 
fC96O2(i)=fC96O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC17H28O5(i)=fC17H28O5(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + C810O2 = C17H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C810O2'; 
fC96O2(i)=fC96O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC17H28O5(i)=fC17H28O5(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + C812O2 = C17H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C812O2'; 
fC96O2(i)=fC96O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC17H28O6(i)=fC17H28O6(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + C813O2 = C17H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'C813O2'; 
fC96O2(i)=fC96O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC17H28O7(i)=fC17H28O7(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC96O2(i)=fC96O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C96O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'CHEXO2'; 
fC96O2(i)=fC96O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C920O2 rxns
i=i+1;
Rnames{i} = 'C920O2 + C920O2 = C18H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'C920O2'; 
fC920O2(i)=fC920O2(i)-1; fC920O2(i)=fC920O2(i)-1; fC18H30O6(i)=fC18H30O6(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + C97O2 = C18H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'C97O2'; 
fC920O2(i)=fC920O2(i)-1; fC97O2(i)=fC97O2(i)-1; fC18H30O6(i)=fC18H30O6(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + C921O2 = C18H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'C921O2'; 
fC920O2(i)=fC920O2(i)-1; fC921O2(i)=fC921O2(i)-1; fC18H30O7(i)=fC18H30O7(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + C98O2 = C18H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'C98O2'; 
fC920O2(i)=fC920O2(i)-1; fC98O2(i)=fC98O2(i)-1; fC18H30O7(i)=fC18H30O7(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + C922O2 = C18H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'C922O2'; 
fC920O2(i)=fC920O2(i)-1; fC922O2(i)=fC922O2(i)-1; fC18H30O8(i)=fC18H30O8(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + C721CO3 = C17H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'C721CO3'; 
fC920O2(i)=fC920O2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC17H26O7(i)=fC17H26O7(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + C85O2 = C17H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'C85O2'; 
fC920O2(i)=fC920O2(i)-1; fC85O2(i)=fC85O2(i)-1; fC17H28O5(i)=fC17H28O5(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + C89O2 = C17H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'C89O2'; 
fC920O2(i)=fC920O2(i)-1; fC89O2(i)=fC89O2(i)-1; fC17H28O5(i)=fC17H28O5(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + C86O2 = C17H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'C86O2'; 
fC920O2(i)=fC920O2(i)-1; fC86O2(i)=fC86O2(i)-1; fC17H28O6(i)=fC17H28O6(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + C811O2 = C17H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'C811O2'; 
fC920O2(i)=fC920O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC17H28O6(i)=fC17H28O6(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + C810O2 = C17H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'C810O2'; 
fC920O2(i)=fC920O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC17H28O6(i)=fC17H28O6(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + C812O2 = C17H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'C812O2'; 
fC920O2(i)=fC920O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC17H28O7(i)=fC17H28O7(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + C813O2 = C17H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'C813O2'; 
fC920O2(i)=fC920O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC17H28O8(i)=fC17H28O8(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC920O2(i)=fC920O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C920O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'CHEXO2'; 
fC920O2(i)=fC920O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C97O2 rxns
i=i+1;
Rnames{i} = 'C97O2 + C97O2 = C18H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'C97O2'; 
fC97O2(i)=fC97O2(i)-1; fC97O2(i)=fC97O2(i)-1; fC18H30O6(i)=fC18H30O6(i)+1;

i=i+1;
Rnames{i} = 'C97O2 + C921O2 = C18H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'C921O2'; 
fC97O2(i)=fC97O2(i)-1; fC921O2(i)=fC921O2(i)-1; fC18H30O7(i)=fC18H30O7(i)+1;

i=i+1;
Rnames{i} = 'C97O2 + C98O2 = C18H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'C98O2'; 
fC97O2(i)=fC97O2(i)-1; fC98O2(i)=fC98O2(i)-1; fC18H30O7(i)=fC18H30O7(i)+1;

i=i+1;
Rnames{i} = 'C97O2 + C922O2 = C18H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'C922O2'; 
fC97O2(i)=fC97O2(i)-1; fC922O2(i)=fC922O2(i)-1; fC18H30O8(i)=fC18H30O8(i)+1;

i=i+1;
Rnames{i} = 'C97O2 + C721CO3 = C17H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'C721CO3'; 
fC97O2(i)=fC97O2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC17H26O7(i)=fC17H26O7(i)+1;

i=i+1;
Rnames{i} = 'C97O2 + C85O2 = C17H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'C85O2'; 
fC97O2(i)=fC97O2(i)-1; fC85O2(i)=fC85O2(i)-1; fC17H28O5(i)=fC17H28O5(i)+1;

i=i+1;
Rnames{i} = 'C97O2 + C89O2 = C17H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'C89O2'; 
fC97O2(i)=fC97O2(i)-1; fC89O2(i)=fC89O2(i)-1; fC17H28O5(i)=fC17H28O5(i)+1;

i=i+1;
Rnames{i} = 'C97O2 + C86O2 = C17H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'C86O2'; 
fC97O2(i)=fC97O2(i)-1; fC86O2(i)=fC86O2(i)-1; fC17H28O6(i)=fC17H28O6(i)+1;

i=i+1;
Rnames{i} = 'C97O2 + C811O2 = C17H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'C811O2'; 
fC97O2(i)=fC97O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC17H28O6(i)=fC17H28O6(i)+1;

i=i+1;
Rnames{i} = 'C97O2 + C810O2 = C17H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'C810O2'; 
fC97O2(i)=fC97O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC17H28O6(i)=fC17H28O6(i)+1;

i=i+1;
Rnames{i} = 'C97O2 + C812O2 = C17H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'C812O2'; 
fC97O2(i)=fC97O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC17H28O7(i)=fC17H28O7(i)+1;

i=i+1;
Rnames{i} = 'C97O2 + C813O2 = C17H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'C813O2'; 
fC97O2(i)=fC97O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC17H28O8(i)=fC17H28O8(i)+1;

i=i+1;
Rnames{i} = 'C97O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC97O2(i)=fC97O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C97O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'CHEXO2'; 
fC97O2(i)=fC97O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C921O2 rxns
i=i+1;
Rnames{i} = 'C921O2 + C921O2 = C18H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'C921O2'; 
fC921O2(i)=fC921O2(i)-1; fC921O2(i)=fC921O2(i)-1; fC18H30O8(i)=fC18H30O8(i)+1;

i=i+1;
Rnames{i} = 'C921O2 + C98O2 = C18H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'C98O2'; 
fC921O2(i)=fC921O2(i)-1; fC98O2(i)=fC98O2(i)-1; fC18H30O8(i)=fC18H30O8(i)+1;

i=i+1;
Rnames{i} = 'C921O2 + C922O2 = C18H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'C922O2'; 
fC921O2(i)=fC921O2(i)-1; fC922O2(i)=fC922O2(i)-1; fC18H30O9(i)=fC18H30O9(i)+1;

i=i+1;
Rnames{i} = 'C921O2 + C721CO3 = C17H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'C721CO3'; 
fC921O2(i)=fC921O2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC17H26O8(i)=fC17H26O8(i)+1;

i=i+1;
Rnames{i} = 'C921O2 + C85O2 = C17H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'C85O2'; 
fC921O2(i)=fC921O2(i)-1; fC85O2(i)=fC85O2(i)-1; fC17H28O6(i)=fC17H28O6(i)+1;

i=i+1;
Rnames{i} = 'C921O2 + C89O2 = C17H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'C89O2'; 
fC921O2(i)=fC921O2(i)-1; fC89O2(i)=fC89O2(i)-1; fC17H28O6(i)=fC17H28O6(i)+1;

i=i+1;
Rnames{i} = 'C921O2 + C86O2 = C17H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'C86O2'; 
fC921O2(i)=fC921O2(i)-1; fC86O2(i)=fC86O2(i)-1; fC17H28O7(i)=fC17H28O7(i)+1;

i=i+1;
Rnames{i} = 'C921O2 + C811O2 = C17H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'C811O2'; 
fC921O2(i)=fC921O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC17H28O7(i)=fC17H28O7(i)+1;

i=i+1;
Rnames{i} = 'C921O2 + C810O2 = C17H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'C810O2'; 
fC921O2(i)=fC921O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC17H28O7(i)=fC17H28O7(i)+1;

i=i+1;
Rnames{i} = 'C921O2 + C812O2 = C17H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'C812O2'; 
fC921O2(i)=fC921O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC17H28O8(i)=fC17H28O8(i)+1;

i=i+1;
Rnames{i} = 'C921O2 + C813O2 = C17H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'C813O2'; 
fC921O2(i)=fC921O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC17H28O9(i)=fC17H28O9(i)+1;

i=i+1;
Rnames{i} = 'C921O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC921O2(i)=fC921O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C921O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'CHEXO2'; 
fC921O2(i)=fC921O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C98O2 rxns
i=i+1;
Rnames{i} = 'C98O2 + C98O2 = C18H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'C98O2'; 
fC98O2(i)=fC98O2(i)-1; fC98O2(i)=fC98O2(i)-1; fC18H30O8(i)=fC18H30O8(i)+1;

i=i+1;
Rnames{i} = 'C98O2 + C922O2 = C18H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'C922O2'; 
fC98O2(i)=fC98O2(i)-1; fC922O2(i)=fC922O2(i)-1; fC18H30O9(i)=fC18H30O9(i)+1;

i=i+1;
Rnames{i} = 'C98O2 + C721CO3 = C17H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'C721CO3'; 
fC98O2(i)=fC98O2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC17H26O8(i)=fC17H26O8(i)+1;

i=i+1;
Rnames{i} = 'C98O2 + C85O2 = C17H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'C85O2'; 
fC98O2(i)=fC98O2(i)-1; fC85O2(i)=fC85O2(i)-1; fC17H28O6(i)=fC17H28O6(i)+1;

i=i+1;
Rnames{i} = 'C98O2 + C89O2 = C17H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'C89O2'; 
fC98O2(i)=fC98O2(i)-1; fC89O2(i)=fC89O2(i)-1; fC17H28O6(i)=fC17H28O6(i)+1;

i=i+1;
Rnames{i} = 'C98O2 + C86O2 = C17H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'C86O2'; 
fC98O2(i)=fC98O2(i)-1; fC86O2(i)=fC86O2(i)-1; fC17H28O7(i)=fC17H28O7(i)+1;

i=i+1;
Rnames{i} = 'C98O2 + C811O2 = C17H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'C811O2'; 
fC98O2(i)=fC98O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC17H28O7(i)=fC17H28O7(i)+1;

i=i+1;
Rnames{i} = 'C98O2 + C810O2 = C17H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'C810O2'; 
fC98O2(i)=fC98O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC17H28O7(i)=fC17H28O7(i)+1;

i=i+1;
Rnames{i} = 'C98O2 + C812O2 = C17H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'C812O2'; 
fC98O2(i)=fC98O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC17H28O8(i)=fC17H28O8(i)+1;

i=i+1;
Rnames{i} = 'C98O2 + C813O2 = C17H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'C813O2'; 
fC98O2(i)=fC98O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC17H28O9(i)=fC17H28O9(i)+1;

i=i+1;
Rnames{i} = 'C98O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC98O2(i)=fC98O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C98O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'CHEXO2'; 
fC98O2(i)=fC98O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C922O2 rxns
i=i+1;
Rnames{i} = 'C922O2 + C922O2 = C18H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'C922O2'; 
fC922O2(i)=fC922O2(i)-1; fC922O2(i)=fC922O2(i)-1; fC18H30O10(i)=fC18H30O10(i)+1;

i=i+1;
Rnames{i} = 'C922O2 + C721CO3 = C17H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'C721CO3'; 
fC922O2(i)=fC922O2(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC17H26O9(i)=fC17H26O9(i)+1;

i=i+1;
Rnames{i} = 'C922O2 + C85O2 = C17H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'C85O2'; 
fC922O2(i)=fC922O2(i)-1; fC85O2(i)=fC85O2(i)-1; fC17H28O7(i)=fC17H28O7(i)+1;

i=i+1;
Rnames{i} = 'C922O2 + C89O2 = C17H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'C89O2'; 
fC922O2(i)=fC922O2(i)-1; fC89O2(i)=fC89O2(i)-1; fC17H28O7(i)=fC17H28O7(i)+1;

i=i+1;
Rnames{i} = 'C922O2 + C86O2 = C17H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'C86O2'; 
fC922O2(i)=fC922O2(i)-1; fC86O2(i)=fC86O2(i)-1; fC17H28O8(i)=fC17H28O8(i)+1;

i=i+1;
Rnames{i} = 'C922O2 + C811O2 = C17H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'C811O2'; 
fC922O2(i)=fC922O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC17H28O8(i)=fC17H28O8(i)+1;

i=i+1;
Rnames{i} = 'C922O2 + C810O2 = C17H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'C810O2'; 
fC922O2(i)=fC922O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC17H28O8(i)=fC17H28O8(i)+1;

i=i+1;
Rnames{i} = 'C922O2 + C812O2 = C17H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'C812O2'; 
fC922O2(i)=fC922O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC17H28O9(i)=fC17H28O9(i)+1;

i=i+1;
Rnames{i} = 'C922O2 + C813O2 = C17H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'C813O2'; 
fC922O2(i)=fC922O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC17H28O10(i)=fC17H28O10(i)+1;

i=i+1;
Rnames{i} = 'C922O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORb.*fROOR1;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC922O2(i)=fC922O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C922O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'CHEXO2'; 
fC922O2(i)=fC922O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C721CO3 rxns
i=i+1;
Rnames{i} = 'C721CO3 + C721CO3 = C16H22O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'C721CO3'; 
fC721CO3(i)=fC721CO3(i)-1; fC721CO3(i)=fC721CO3(i)-1; fC16H22O8(i)=fC16H22O8(i)+1;

i=i+1;
Rnames{i} = 'C721CO3 + C85O2 = C16H24O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'C85O2'; 
fC721CO3(i)=fC721CO3(i)-1; fC85O2(i)=fC85O2(i)-1; fC16H24O6(i)=fC16H24O6(i)+1;

i=i+1;
Rnames{i} = 'C721CO3 + C89O2 = C16H24O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'C89O2'; 
fC721CO3(i)=fC721CO3(i)-1; fC89O2(i)=fC89O2(i)-1; fC16H24O6(i)=fC16H24O6(i)+1;

i=i+1;
Rnames{i} = 'C721CO3 + C86O2 = C16H24O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'C86O2'; 
fC721CO3(i)=fC721CO3(i)-1; fC86O2(i)=fC86O2(i)-1; fC16H24O7(i)=fC16H24O7(i)+1;

i=i+1;
Rnames{i} = 'C721CO3 + C811O2 = C16H24O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'C811O2'; 
fC721CO3(i)=fC721CO3(i)-1; fC811O2(i)=fC811O2(i)-1; fC16H24O7(i)=fC16H24O7(i)+1;

i=i+1;
Rnames{i} = 'C721CO3 + C810O2 = C16H24O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'C810O2'; 
fC721CO3(i)=fC721CO3(i)-1; fC810O2(i)=fC810O2(i)-1; fC16H24O7(i)=fC16H24O7(i)+1;

i=i+1;
Rnames{i} = 'C721CO3 + C812O2 = C16H24O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'C812O2'; 
fC721CO3(i)=fC721CO3(i)-1; fC812O2(i)=fC812O2(i)-1; fC16H24O8(i)=fC16H24O8(i)+1;

i=i+1;
Rnames{i} = 'C721CO3 + C813O2 = C16H24O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'C813O2'; 
fC721CO3(i)=fC721CO3(i)-1; fC813O2(i)=fC813O2(i)-1; fC16H24O9(i)=fC16H24O9(i)+1;

i=i+1;
Rnames{i} = 'C721CO3 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'BUT2OLO2'; 
fC721CO3(i)=fC721CO3(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C721CO3 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'CHEXO2'; 
fC721CO3(i)=fC721CO3(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C85O2 rxns
i=i+1;
Rnames{i} = 'C85O2 + C85O2 = C16H26O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'C85O2'; 
fC85O2(i)=fC85O2(i)-1; fC85O2(i)=fC85O2(i)-1; fC16H26O4(i)=fC16H26O4(i)+1;

i=i+1;
Rnames{i} = 'C85O2 + C89O2 = C16H26O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'C89O2'; 
fC85O2(i)=fC85O2(i)-1; fC89O2(i)=fC89O2(i)-1; fC16H26O4(i)=fC16H26O4(i)+1;

i=i+1;
Rnames{i} = 'C85O2 + C86O2 = C16H26O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'C86O2'; 
fC85O2(i)=fC85O2(i)-1; fC86O2(i)=fC86O2(i)-1; fC16H26O5(i)=fC16H26O5(i)+1;

i=i+1;
Rnames{i} = 'C85O2 + C811O2 = C16H26O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'C811O2'; 
fC85O2(i)=fC85O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC16H26O5(i)=fC16H26O5(i)+1;

i=i+1;
Rnames{i} = 'C85O2 + C810O2 = C16H26O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'C810O2'; 
fC85O2(i)=fC85O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC16H26O5(i)=fC16H26O5(i)+1;

i=i+1;
Rnames{i} = 'C85O2 + C812O2 = C16H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'C812O2'; 
fC85O2(i)=fC85O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC16H26O6(i)=fC16H26O6(i)+1;

i=i+1;
Rnames{i} = 'C85O2 + C813O2 = C16H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'C813O2'; 
fC85O2(i)=fC85O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC16H26O7(i)=fC16H26O7(i)+1;

i=i+1;
Rnames{i} = 'C85O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC85O2(i)=fC85O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C85O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'CHEXO2'; 
fC85O2(i)=fC85O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C89O2 rxns
i=i+1;
Rnames{i} = 'C89O2 + C89O2 = C16H26O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'C89O2'; 
fC89O2(i)=fC89O2(i)-1; fC89O2(i)=fC89O2(i)-1; fC16H26O4(i)=fC16H26O4(i)+1;

i=i+1;
Rnames{i} = 'C89O2 + C86O2 = C16H26O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'C86O2'; 
fC89O2(i)=fC89O2(i)-1; fC86O2(i)=fC86O2(i)-1; fC16H26O5(i)=fC16H26O5(i)+1;

i=i+1;
Rnames{i} = 'C89O2 + C811O2 = C16H26O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'C811O2'; 
fC89O2(i)=fC89O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC16H26O5(i)=fC16H26O5(i)+1;

i=i+1;
Rnames{i} = 'C89O2 + C810O2 = C16H26O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'C810O2'; 
fC89O2(i)=fC89O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC16H26O5(i)=fC16H26O5(i)+1;

i=i+1;
Rnames{i} = 'C89O2 + C812O2 = C16H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'C812O2'; 
fC89O2(i)=fC89O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC16H26O6(i)=fC16H26O6(i)+1;

i=i+1;
Rnames{i} = 'C89O2 + C813O2 = C16H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'C813O2'; 
fC89O2(i)=fC89O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC16H26O7(i)=fC16H26O7(i)+1;

i=i+1;
Rnames{i} = 'C89O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC89O2(i)=fC89O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C89O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'CHEXO2'; 
fC89O2(i)=fC89O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C86O2 rxns
i=i+1;
Rnames{i} = 'C86O2 + C86O2 = C16H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'C86O2'; 
fC86O2(i)=fC86O2(i)-1; fC86O2(i)=fC86O2(i)-1; fC16H26O6(i)=fC16H26O6(i)+1;

i=i+1;
Rnames{i} = 'C86O2 + C811O2 = C16H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'C811O2'; 
fC86O2(i)=fC86O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC16H26O6(i)=fC16H26O6(i)+1;

i=i+1;
Rnames{i} = 'C86O2 + C810O2 = C16H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'C810O2'; 
fC86O2(i)=fC86O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC16H26O6(i)=fC16H26O6(i)+1;

i=i+1;
Rnames{i} = 'C86O2 + C812O2 = C16H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'C812O2'; 
fC86O2(i)=fC86O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC16H26O7(i)=fC16H26O7(i)+1;

i=i+1;
Rnames{i} = 'C86O2 + C813O2 = C16H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'C813O2'; 
fC86O2(i)=fC86O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC16H26O8(i)=fC16H26O8(i)+1;

i=i+1;
Rnames{i} = 'C86O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC86O2(i)=fC86O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C86O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'CHEXO2'; 
fC86O2(i)=fC86O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C811O2 rxns
i=i+1;
Rnames{i} = 'C811O2 + C811O2 = C16H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'C811O2'; 
fC811O2(i)=fC811O2(i)-1; fC811O2(i)=fC811O2(i)-1; fC16H26O6(i)=fC16H26O6(i)+1;

i=i+1;
Rnames{i} = 'C811O2 + C810O2 = C16H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'C810O2'; 
fC811O2(i)=fC811O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC16H26O6(i)=fC16H26O6(i)+1;

i=i+1;
Rnames{i} = 'C811O2 + C812O2 = C16H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'C812O2'; 
fC811O2(i)=fC811O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC16H26O7(i)=fC16H26O7(i)+1;

i=i+1;
Rnames{i} = 'C811O2 + C813O2 = C16H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'C813O2'; 
fC811O2(i)=fC811O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC16H26O8(i)=fC16H26O8(i)+1;

i=i+1;
Rnames{i} = 'C811O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC811O2(i)=fC811O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C811O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'CHEXO2'; 
fC811O2(i)=fC811O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C810O2 rxns
i=i+1;
Rnames{i} = 'C810O2 + C810O2 = C16H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'C810O2'; 
fC810O2(i)=fC810O2(i)-1; fC810O2(i)=fC810O2(i)-1; fC16H26O6(i)=fC16H26O6(i)+1;

i=i+1;
Rnames{i} = 'C810O2 + C812O2 = C16H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'C812O2'; 
fC810O2(i)=fC810O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC16H26O7(i)=fC16H26O7(i)+1;

i=i+1;
Rnames{i} = 'C810O2 + C813O2 = C16H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'C813O2'; 
fC810O2(i)=fC810O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC16H26O8(i)=fC16H26O8(i)+1;

i=i+1;
Rnames{i} = 'C810O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC810O2(i)=fC810O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C810O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'CHEXO2'; 
fC810O2(i)=fC810O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C812O2 rxns
i=i+1;
Rnames{i} = 'C812O2 + C812O2 = C16H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'C812O2'; 
fC812O2(i)=fC812O2(i)-1; fC812O2(i)=fC812O2(i)-1; fC16H26O8(i)=fC16H26O8(i)+1;

i=i+1;
Rnames{i} = 'C812O2 + C813O2 = C16H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'C813O2'; 
fC812O2(i)=fC812O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC16H26O9(i)=fC16H26O9(i)+1;

i=i+1;
Rnames{i} = 'C812O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC812O2(i)=fC812O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C812O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'CHEXO2'; 
fC812O2(i)=fC812O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

% C813O2 rxns
i=i+1;
Rnames{i} = 'C813O2 + C813O2 = C16H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'C813O2'; 
fC813O2(i)=fC813O2(i)-1; fC813O2(i)=fC813O2(i)-1; fC16H26O10(i)=fC16H26O10(i)+1;

i=i+1;
Rnames{i} = 'C813O2 + BUT2OLO2 = BUT2OLdimer';
k(:,i) = kROORb.*fROOR1;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'BUT2OLO2'; 
fC813O2(i)=fC813O2(i)-1; fBUT2OLO2(i)=fBUT2OLO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1;

i=i+1;
Rnames{i} = 'C813O2 + CHEXO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'CHEXO2'; 
fC813O2(i)=fC813O2(i)-1; fCHEXO2(i)=fCHEXO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1;

%% Dimers formation between new APRO2s (non-NO3)
% PINAL4RO2 rxns
i=i+1;
Rnames{i} = 'PINAL4RO2 + PINAL4RO2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'PINAL4RO2'; 
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1;

i=i+1;
Rnames{i} = 'PINAL4RO2 + APC10H15O4RB = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H30O6(i)=fC20H30O6(i)+1;

i=i+1;
Rnames{i} = 'PINAL4RO2 + PINAL4H1RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'PINAL4RO2 + APC10H15O5RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'PINAL4RO2 + PINALPA4RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'PINAL4RO2 + APC10H15O6RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'PINAL4RO2 + PINAL10P1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'PINAL4RO2 + PINALPA1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'PINAL4RO2 + APC10H15O7RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1;

i=i+1;
Rnames{i} = 'PINAL4RO2 + APC10H15O8RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'PINAL4RO2 + APC10H15O9RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1;

i=i+1;
Rnames{i} = 'PINAL4RO2 + APC10H15O10RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINAL4RO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fPINAL4RO2(i)=fPINAL4RO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1;

% APC10H15O4RB rxns
i=i+1;
Rnames{i} = 'APC10H15O4RB + APC10H15O4RB = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'APC10H15O4RB'; 
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H30O6(i)=fC20H30O6(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O4RB + PINAL4H1RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'PINAL4H1RO2'; 
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O4RB + APC10H15O5RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'APC10H15O5RO2'; 
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O4RB + PINALPA4RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'PINALPA4RO2'; 
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O4RB + APC10H15O6RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'APC10H15O6RO2'; 
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O4RB + PINAL10P1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'PINAL10P1RO2'; 
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O4RB + PINALPA1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'PINALPA1RO2'; 
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O4RB + APC10H15O7RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'APC10H15O7RO2'; 
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O4RB + APC10H15O8RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'APC10H15O8RO2'; 
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O4RB + APC10H15O9RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'APC10H15O9RO2'; 
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O4RB + APC10H15O10RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APC10H15O4RB'; Gstr{i,2} = 'APC10H15O10RO2'; 
fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1;


% PINAL4H1RO2 rxns
i=i+1;
Rnames{i} = 'PINAL4H1RO2 + PINAL4H1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'PINAL4H1RO2 + APC10H15O5RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'PINAL4H1RO2 + PINALPA4RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1;

i=i+1;
Rnames{i} = 'PINAL4H1RO2 + APC10H15O6RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1;

i=i+1;
Rnames{i} = 'PINAL4H1RO2 + PINAL10P1RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1;

i=i+1;
Rnames{i} = 'PINAL4H1RO2 + PINALPA1RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1;

i=i+1;
Rnames{i} = 'PINAL4H1RO2 + APC10H15O7RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'PINAL4H1RO2 + APC10H15O8RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1;

i=i+1;
Rnames{i} = 'PINAL4H1RO2 + APC10H15O9RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1;

i=i+1;
Rnames{i} = 'PINAL4H1RO2 + APC10H15O10RO2 = C20H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'PINAL4H1RO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O13(i)=fC20H30O13(i)+1;


% APC10H15O5RO2 rxns
i=i+1;
Rnames{i} = 'APC10H15O5RO2 + APC10H15O5RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O5RO2 + PINALPA4RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O5RO2 + APC10H15O6RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O5RO2 + PINAL10P1RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O5RO2 + PINALPA1RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O5RO2 + APC10H15O7RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O5RO2 + APC10H15O8RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O5RO2 + APC10H15O9RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O5RO2 + APC10H15O10RO2 = C20H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'APC10H15O5RO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O13(i)=fC20H30O13(i)+1;

% PINALPA4RO2 rxns
i=i+1;
Rnames{i} = 'PINALPA4RO2 + PINALPA4RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'PINALPA4RO2 + APC10H15O6RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'PINALPA4RO2 + PINAL10P1RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'PINALPA4RO2 + PINALPA1RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'PINALPA4RO2 + APC10H15O7RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1;

i=i+1;
Rnames{i} = 'PINALPA4RO2 + APC10H15O8RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1;

i=i+1;
Rnames{i} = 'PINALPA4RO2 + APC10H15O9RO2 = C20H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O13(i)=fC20H30O13(i)+1;

i=i+1;
Rnames{i} = 'PINALPA4RO2 + APC10H15O10RO2 = C20H30O14';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'PINALPA4RO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O14(i)=fC20H30O14(i)+1;

% APC10H15O6RO2 rxns
i=i+1;
Rnames{i} = 'APC10H15O6RO2 + APC10H15O6RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O6RO2 + PINAL10P1RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O6RO2 + PINALPA1RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O6RO2 + APC10H15O7RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O6RO2 + APC10H15O8RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O6RO2 + APC10H15O9RO2 = C20H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O13(i)=fC20H30O13(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O6RO2 + APC10H15O10RO2 = C20H30O14';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'APC10H15O6RO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O14(i)=fC20H30O14(i)+1;

% PINAL10P1RO2 rxns
i=i+1;
Rnames{i} = 'PINAL10P1RO2 + PINAL10P1RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINAL10P1RO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'PINAL10P1RO2 + PINALPA1RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINAL10P1RO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'PINAL10P1RO2 + APC10H15O7RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINAL10P1RO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1;

i=i+1;
Rnames{i} = 'PINAL10P1RO2 + APC10H15O8RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINAL10P1RO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1;

i=i+1;
Rnames{i} = 'PINAL10P1RO2 + APC10H15O9RO2 = C20H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'PINAL10P1RO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O13(i)=fC20H30O13(i)+1;

i=i+1;
Rnames{i} = 'PINAL10P1RO2 + APC10H15O10RO2 = C20H30O14';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'PINAL10P1RO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O14(i)=fC20H30O14(i)+1;

% PINALPA1RO2 rxns
i=i+1;
Rnames{i} = 'PINALPA1RO2 + PINALPA1RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINALPA1RO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1;

i=i+1;
Rnames{i} = 'PINALPA1RO2 + APC10H15O7RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINALPA1RO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1;

i=i+1;
Rnames{i} = 'PINALPA1RO2 + APC10H15O8RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINALPA1RO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1;

i=i+1;
Rnames{i} = 'PINALPA1RO2 + APC10H15O9RO2 = C20H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'PINALPA1RO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O13(i)=fC20H30O13(i)+1;

i=i+1;
Rnames{i} = 'PINALPA1RO2 + APC10H15O10RO2 = C20H30O14';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'PINALPA1RO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O14(i)=fC20H30O14(i)+1;

% APC10H15O7RO2 rxns
i=i+1;
Rnames{i} = 'APC10H15O7RO2 + APC10H15O7RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APC10H15O7RO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O7RO2 + APC10H15O8RO2 = C20H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'APC10H15O7RO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O13(i)=fC20H30O13(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O7RO2 + APC10H15O9RO2 = C20H30O14';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'APC10H15O7RO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O14(i)=fC20H30O14(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O7RO2 + APC10H15O10RO2 = C20H30O15';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'APC10H15O7RO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O15(i)=fC20H30O15(i)+1;

% APC10H15O8RO2 rxns
i=i+1;
Rnames{i} = 'APC10H15O8RO2 + APC10H15O8RO2 = C20H30O14';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'APC10H15O8RO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O14(i)=fC20H30O14(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O8RO2 + APC10H15O9RO2 = C20H30O15';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'APC10H15O8RO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O15(i)=fC20H30O15(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O8RO2 + APC10H15O10RO2 = C20H30O16';
k(:,i) = kROORe.*fROOR;
Gstr{i,1} = 'APC10H15O8RO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O16(i)=fC20H30O16(i)+1;

% APC10H15O9RO2 rxns
i=i+1;
Rnames{i} = 'APC10H15O9RO2 + APC10H15O9RO2 = C20H30O16';
k(:,i) = kROORe.*fROOR;
Gstr{i,1} = 'APC10H15O9RO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O16(i)=fC20H30O16(i)+1;

i=i+1;
Rnames{i} = 'APC10H15O9RO2 + APC10H15O10RO2 = C20H30O17';
k(:,i) = kROORe.*fROOR;
Gstr{i,1} = 'APC10H15O9RO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O17(i)=fC20H30O17(i)+1;

% APC10H15O10RO2 rxns
i=i+1;
Rnames{i} = 'APC10H15O10RO2 + APC10H15O10RO2 = C20H30O18';
k(:,i) = kROORe.*fROOR;
Gstr{i,1} = 'APC10H15O10RO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O18(i)=fC20H30O18(i)+1;

%% PINAL4RO2 rxns
i=i+1;
Rnames{i} = 'APINAO2 + PINAL4RO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'PINAL4RO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'APINBO2 + PINAL4RO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'PINAL4RO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'APINCO2 + PINAL4RO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'PINAL4RO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'C107O2 + PINAL4RO2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC107O2(i)=fC107O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C109O2 + PINAL4RO2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC109O2(i)=fC109O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C106O2 + PINAL4RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC106O2(i)=fC106O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C920CO3 + PINAL4RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC920CO3(i)=fC920CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C108O2 + PINAL4RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC108O2(i)=fC108O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'PINALO2 + PINAL4RO2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'PINAL4RO2'; 
fPINALO2(i)=fPINALO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C96CO3 + PINAL4RO2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC96CO3(i)=fC96CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C923CO3 + PINAL4RO2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C923CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC923CO3(i)=fC923CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'LIMAO2 + PINAL4RO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMAO2'; Gstr{i,2} = 'PINAL4RO2'; 
fLIMAO2(i)=fLIMAO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'LIMALBO2 + PINAL4RO2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMALBO2'; Gstr{i,2} = 'PINAL4RO2'; 
fLIMALBO2(i)=fLIMALBO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'LIMCO2 + PINAL4RO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMCO2'; Gstr{i,2} = 'PINAL4RO2'; 
fLIMCO2(i)=fLIMCO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'LIMALO2 + PINAL4RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALO2'; Gstr{i,2} = 'PINAL4RO2'; 
fLIMALO2(i)=fLIMALO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H30O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMBO2 + PINAL4RO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMBO2'; Gstr{i,2} = 'PINAL4RO2'; 
fLIMBO2(i)=fLIMBO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'LIMALAO2 + PINAL4RO2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMALAO2'; Gstr{i,2} = 'PINAL4RO2'; 
fLIMALAO2(i)=fLIMALAO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'NAPINAO2 + PINAL4RO2 = C20H31O4NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINAO2'; Gstr{i,2} = 'PINAL4RO2'; 
fNAPINAO2(i)=fNAPINAO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H31O4NO3(i)=fC20H31O4NO3(i)+1; 

i=i+1;
Rnames{i} = 'NAPINBO2 + PINAL4RO2 = C20H31O4NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINBO2'; Gstr{i,2} = 'PINAL4RO2'; 
fNAPINBO2(i)=fNAPINBO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H31O4NO3(i)=fC20H31O4NO3(i)+1; 

i=i+1;
Rnames{i} = 'BPINAO2 + PINAL4RO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'BPINAO2'; Gstr{i,2} = 'PINAL4RO2'; 
fBPINAO2(i)=fBPINAO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'BPINBO2 + PINAL4RO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'BPINBO2'; Gstr{i,2} = 'PINAL4RO2'; 
fBPINBO2(i)=fBPINBO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H32O5(i)=fC20H30O5(i)+1; 

i=i+1;
Rnames{i} = 'BPINCO2 + PINAL4RO2 = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'BPINCO2'; Gstr{i,2} = 'PINAL4RO2'; 
fBPINCO2(i)=fBPINCO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'C918CO3 + PINAL4RO2 = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C918CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC918CO3(i)=fC918CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'NC102O2 + PINAL4RO2 = C20H29O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NC102O2'; Gstr{i,2} = 'PINAL4RO2'; 
fNC102O2(i)=fNC102O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H29O6NO3(i)=fC20H29O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC101O2 + PINAL4RO2 = C20H29O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NC101O2'; Gstr{i,2} = 'PINAL4RO2'; 
fNC101O2(i)=fNC101O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H29O5NO3(i)=fC20H29O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMO2 + PINAL4RO2 = C20H31O4NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLIMO2'; Gstr{i,2} = 'PINAL4RO2'; 
fNLIMO2(i)=fNLIMO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H31O4NO3(i)=fC20H31O4NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMALO2 + PINAL4RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLIMALO2'; Gstr{i,2} = 'PINAL4RO2'; 
fNLIMALO2(i)=fNLIMALO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC91CO3 + PINAL4RO2 = C20H29O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NC91CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fNC91CO3(i)=fNC91CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H29O5NO3(i)=fC20H29O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINAO2 + PINAL4RO2 = C20H31O4NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINAO2'; Gstr{i,2} = 'PINAL4RO2'; 
fNBPINAO2(i)=fNBPINAO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H31O4NO3(i)=fC20H31O4NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINBO2 + PINAL4RO2 = C20H31O4NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINBO2'; Gstr{i,2} = 'PINAL4RO2'; 
fNBPINBO2(i)=fNBPINBO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC20H31O4NO3(i)=fC20H31O4NO3(i)+1; 

i=i+1;
Rnames{i} = 'C96O2 + PINAL4RO2 = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC96O2(i)=fC96O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H30O5(i)=fC19H30O5(i)+1; 

i=i+1;
Rnames{i} = 'C89CO3 + PINAL4RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC89CO3(i)=fC89CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C920O2 + PINAL4RO2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC920O2(i)=fC920O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C97O2 + PINAL4RO2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC97O2(i)=fC97O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C85CO3 + PINAL4RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC85CO3(i)=fC85CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C811CO3 + PINAL4RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC811CO3(i)=fC811CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C921O2 + PINAL4RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC921O2(i)=fC921O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C98O2 + PINAL4RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC98O2(i)=fC98O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C922O2 + PINAL4RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC922O2(i)=fC922O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C923O2 + PINAL4RO2 = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C923O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC923O2(i)=fC923O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H30O5(i)=fC19H30O5(i)+1; 

i=i+1;
Rnames{i} = 'C924O2 + PINAL4RO2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C924O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC924O2(i)=fC924O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C816CO3 + PINAL4RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C816CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC816CO3(i)=fC816CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'NORLIMO2 + PINAL4RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NORLIMO2'; Gstr{i,2} = 'PINAL4RO2'; 
fNORLIMO2(i)=fNORLIMO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'LMKAO2 + PINAL4RO2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LMKAO2'; Gstr{i,2} = 'PINAL4RO2'; 
fLMKAO2(i)=fLMKAO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'LMKBO2 + PINAL4RO2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LMKBO2'; Gstr{i,2} = 'PINAL4RO2'; 
fLMKBO2(i)=fLMKBO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C926O2 + PINAL4RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C926O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC926O2(i)=fC926O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C817CO3 + PINAL4RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC817CO3(i)=fC817CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'LMLKAO2 + PINAL4RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKAO2'; Gstr{i,2} = 'PINAL4RO2'; 
fLMLKAO2(i)=fLMLKAO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'LMLKBO2 + PINAL4RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKBO2'; Gstr{i,2} = 'PINAL4RO2'; 
fLMLKBO2(i)=fLMLKBO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C823CO3 + PINAL4RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC823CO3(i)=fC823CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C925O2 + PINAL4RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C925O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC925O2(i)=fC925O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'NOPINAO2 + PINAL4RO2 = C19H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINAO2'; Gstr{i,2} = 'PINAL4RO2'; 
fNOPINAO2(i)=fNOPINAO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O5(i)=fC19H28O5(i)+1; 

i=i+1;
Rnames{i} = 'NOPINBO2 + PINAL4RO2 = C19H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINBO2'; Gstr{i,2} = 'PINAL4RO2'; 
fNOPINBO2(i)=fNOPINBO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O5(i)=fC19H28O5(i)+1; 

i=i+1;
Rnames{i} = 'NOPINCO2 + PINAL4RO2 = C19H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINCO2'; Gstr{i,2} = 'PINAL4RO2'; 
fNOPINCO2(i)=fNOPINCO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O5(i)=fC19H28O5(i)+1; 

i=i+1;
Rnames{i} = 'NOPINDO2 + PINAL4RO2 = C19H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINDO2'; Gstr{i,2} = 'PINAL4RO2'; 
fNOPINDO2(i)=fNOPINDO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O5(i)=fC19H28O5(i)+1; 

i=i+1;
Rnames{i} = 'C918O2 + PINAL4RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C918O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC918O2(i)=fC918O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C9DCO2 + PINAL4RO2 = C19H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C9DCO2'; Gstr{i,2} = 'PINAL4RO2'; 
fC9DCO2(i)=fC9DCO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H26O6(i)=fC19H26O6(i)+1; 

i=i+1;
Rnames{i} = 'C915O2 + PINAL4RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C915O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC915O2(i)=fC915O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C917O2 + PINAL4RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C917O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC917O2(i)=fC917O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C919O2 + PINAL4RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C919O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC919O2(i)=fC919O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C914O2 + PINAL4RO2 = C19H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C914O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC914O2(i)=fC914O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H26O7(i)=fC19H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C916O2 + PINAL4RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C916O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC916O2(i)=fC916O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C88CO3 + PINAL4RO2 = C19H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC88CO3(i)=fC88CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H26O7(i)=fC19H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C87CO3 + PINAL4RO2 = C19H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C87CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC87CO3(i)=fC87CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H26O8(i)=fC19H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C822CO3 + PINAL4RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C822CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC822CO3(i)=fC822CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'NLMKAO2 + PINAL4RO2 = C19H29O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLMKAO2'; Gstr{i,2} = 'PINAL4RO2'; 
fNLMKAO2(i)=fNLMKAO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC19H29O5NO3(i)=fC19H29O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'C85O2 + PINAL4RO2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC85O2(i)=fC85O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1; 

i=i+1;
Rnames{i} = 'C89O2 + PINAL4RO2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC89O2(i)=fC89O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1; 

i=i+1;
Rnames{i} = 'C86O2 + PINAL4RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC86O2(i)=fC86O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C811O2 + PINAL4RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC811O2(i)=fC811O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C810O2 + PINAL4RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC810O2(i)=fC810O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C812O2 + PINAL4RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC812O2(i)=fC812O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C813O2 + PINAL4RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC813O2(i)=fC813O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C729CO3 + PINAL4RO2 = C18H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C729CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC729CO3(i)=fC729CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H26O6(i)=fC18H26O6(i)+1; 

i=i+1;
Rnames{i} = 'C816O2 + PINAL4RO2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C816O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC816O2(i)=fC816O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1; 

i=i+1;
Rnames{i} = 'C817O2 + PINAL4RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C817O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC817O2(i)=fC817O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C826O2 + PINAL4RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C826O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC826O2(i)=fC826O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C822O2 + PINAL4RO2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C822O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC822O2(i)=fC822O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1; 

i=i+1;
Rnames{i} = 'C818O2 + PINAL4RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C818O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC818O2(i)=fC818O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C823O2 + PINAL4RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C823O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC823O2(i)=fC823O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C819O2 + PINAL4RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C819O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC819O2(i)=fC819O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C727CO3 + PINAL4RO2 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C727CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC727CO3(i)=fC727CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C731CO3 + PINAL4RO2 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C731CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC731CO3(i)=fC731CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C824O2 + PINAL4RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C824O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC824O2(i)=fC824O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C820O2 + PINAL4RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C820O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC820O2(i)=fC820O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C825O2 + PINAL4RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C825O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC825O2(i)=fC825O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C732CO3 + PINAL4RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C732CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC732CO3(i)=fC732CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C8BCO2 + PINAL4RO2 = C18H28O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C8BCO2'; Gstr{i,2} = 'PINAL4RO2'; 
fC8BCO2(i)=fC8BCO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H28O4(i)=fC18H28O4(i)+1; 

i=i+1;
Rnames{i} = 'C88O2 + PINAL4RO2 = C18H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C88O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC88O2(i)=fC88O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H26O6(i)=fC18H26O6(i)+1; 

i=i+1;
Rnames{i} = 'C718CO3 + PINAL4RO2 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C718CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC718CO3(i)=fC718CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C87O2 + PINAL4RO2 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C87O2'; Gstr{i,2} = 'PINAL4RO2'; 
fC87O2(i)=fC87O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C721CO3 + PINAL4RO2 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'PINAL4RO2'; 
fC721CO3(i)=fC721CO3(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'NC826O2 + PINAL4RO2 = C18H27O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NC826O2'; Gstr{i,2} = 'PINAL4RO2'; 
fNC826O2(i)=fNC826O2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fC18H27O6NO3(i)=fC18H27O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'BUT2OLO2 + PINAL4RO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'BUT2OLO2'; Gstr{i,2} = 'PINAL4RO2'; 
fBUT2OLO2(i)=fBUT2OLO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1; 

i=i+1;
Rnames{i} = 'CHEXO2 + PINAL4RO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'CHEXO2'; Gstr{i,2} = 'PINAL4RO2'; 
fCHEXO2(i)=fCHEXO2(i)-1; fPINAL4RO2(i)=fPINAL4RO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1; 

%% APC10H15O4RB rxns
i=i+1;
Rnames{i} = 'APINAO2 + APC10H15O4RB = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fAPINAO2(i)=fAPINAO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'APINBO2 + APC10H15O4RB = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fAPINBO2(i)=fAPINBO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'APINCO2 + APC10H15O4RB = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fAPINCO2(i)=fAPINCO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'C107O2 + APC10H15O4RB = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC107O2(i)=fC107O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C109O2 + APC10H15O4RB = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC109O2(i)=fC109O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C106O2 + APC10H15O4RB = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC106O2(i)=fC106O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C920CO3 + APC10H15O4RB = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC920CO3(i)=fC920CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C108O2 + APC10H15O4RB = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC108O2(i)=fC108O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'PINALO2 + APC10H15O4RB = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fPINALO2(i)=fPINALO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C96CO3 + APC10H15O4RB = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC96CO3(i)=fC96CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C923CO3 + APC10H15O4RB = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C923CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC923CO3(i)=fC923CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'LIMAO2 + APC10H15O4RB = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMAO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fLIMAO2(i)=fLIMAO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'LIMALBO2 + APC10H15O4RB = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMALBO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fLIMALBO2(i)=fLIMALBO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'LIMCO2 + APC10H15O4RB = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMCO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fLIMCO2(i)=fLIMCO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'LIMALO2 + APC10H15O4RB = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fLIMALO2(i)=fLIMALO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMBO2 + APC10H15O4RB = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMBO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fLIMBO2(i)=fLIMBO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'LIMALAO2 + APC10H15O4RB = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMALAO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fLIMALAO2(i)=fLIMALAO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'NAPINAO2 + APC10H15O4RB = C20H31O4NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINAO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNAPINAO2(i)=fNAPINAO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H31O4NO3(i)=fC20H31O4NO3(i)+1; 

i=i+1;
Rnames{i} = 'NAPINBO2 + APC10H15O4RB = C20H31O4NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINBO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNAPINBO2(i)=fNAPINBO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H31O4NO3(i)=fC20H31O4NO3(i)+1; 

i=i+1;
Rnames{i} = 'BPINAO2 + APC10H15O4RB = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'BPINAO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fBPINAO2(i)=fBPINAO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'BPINBO2 + APC10H15O4RB = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'BPINBO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fBPINBO2(i)=fBPINBO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'BPINCO2 + APC10H15O4RB = C20H32O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'BPINCO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fBPINCO2(i)=fBPINCO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H32O5(i)=fC20H32O5(i)+1; 

i=i+1;
Rnames{i} = 'C918CO3 + APC10H15O4RB = C20H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C918CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC918CO3(i)=fC918CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H30O6(i)=fC20H30O6(i)+1; 

i=i+1;
Rnames{i} = 'NC102O2 + APC10H15O4RB = C20H29O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NC102O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNC102O2(i)=fNC102O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H29O6NO3(i)=fC20H29O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC101O2 + APC10H15O4RB = C20H29O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NC101O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNC101O2(i)=fNC101O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H29O5NO3(i)=fC20H29O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMO2 + APC10H15O4RB = C20H31O4NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLIMO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNLIMO2(i)=fNLIMO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H31O4NO3(i)=fC20H31O4NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMALO2 + APC10H15O4RB = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLIMALO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNLIMALO2(i)=fNLIMALO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC91CO3 + APC10H15O4RB = C20H29O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NC91CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fNC91CO3(i)=fNC91CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H29O5NO3(i)=fC20H29O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINAO2 + APC10H15O4RB = C20H31O4NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINAO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNBPINAO2(i)=fNBPINAO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H31O4NO3(i)=fC20H31O4NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINBO2 + APC10H15O4RB = C20H31O4NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINBO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNBPINBO2(i)=fNBPINBO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC20H31O4NO3(i)=fC20H31O4NO3(i)+1; 

i=i+1;
Rnames{i} = 'C96O2 + APC10H15O4RB = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC96O2(i)=fC96O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H30O5(i)=fC19H30O5(i)+1; 

i=i+1;
Rnames{i} = 'C89CO3 + APC10H15O4RB = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC89CO3(i)=fC89CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C920O2 + APC10H15O4RB = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC920O2(i)=fC920O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H30O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C97O2 + APC10H15O4RB = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC97O2(i)=fC97O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H30O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C85CO3 + APC10H15O4RB = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC85CO3(i)=fC85CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C811CO3 + APC10H15O4RB = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC811CO3(i)=fC811CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C921O2 + APC10H15O4RB = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC921O2(i)=fC921O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C98O2 + APC10H15O4RB = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC98O2(i)=fC98O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C922O2 + APC10H15O4RB = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC922O2(i)=fC922O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C923O2 + APC10H15O4RB = C19H30O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C923O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC923O2(i)=fC923O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H30O5(i)=fC19H30O5(i)+1; 

i=i+1;
Rnames{i} = 'C924O2 + APC10H15O4RB = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C924O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC924O2(i)=fC924O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H30O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C816CO3 + APC10H15O4RB = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C816CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC816CO3(i)=fC816CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'NORLIMO2 + APC10H15O4RB = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NORLIMO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNORLIMO2(i)=fNORLIMO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'LMKAO2 + APC10H15O4RB = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LMKAO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fLMKAO2(i)=fLMKAO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H30O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'LMKBO2 + APC10H15O4RB = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LMKBO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fLMKBO2(i)=fLMKBO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H30O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C926O2 + APC10H15O4RB = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C926O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC926O2(i)=fC926O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C817CO3 + APC10H15O4RB = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC817CO3(i)=fC817CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'LMLKAO2 + APC10H15O4RB = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKAO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fLMLKAO2(i)=fLMLKAO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'LMLKBO2 + APC10H15O4RB = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKBO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fLMLKBO2(i)=fLMLKBO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C823CO3 + APC10H15O4RB = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC823CO3(i)=fC823CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C925O2 + APC10H15O4RB = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C925O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC925O2(i)=fC925O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'NOPINAO2 + APC10H15O4RB = C19H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINAO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNOPINAO2(i)=fNOPINAO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O5(i)=fC19H28O5(i)+1; 

i=i+1;
Rnames{i} = 'NOPINBO2 + APC10H15O4RB = C19H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINBO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNOPINBO2(i)=fNOPINBO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O5(i)=fC19H28O5(i)+1; 

i=i+1;
Rnames{i} = 'NOPINCO2 + APC10H15O4RB = C19H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINCO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNOPINCO2(i)=fNOPINCO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O5(i)=fC19H28O5(i)+1; 

i=i+1;
Rnames{i} = 'NOPINDO2 + APC10H15O4RB = C19H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINDO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNOPINDO2(i)=fNOPINDO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O5(i)=fC19H28O5(i)+1; 

i=i+1;
Rnames{i} = 'C918O2 + APC10H15O4RB = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C918O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC918O2(i)=fC918O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C9DCO2 + APC10H15O4RB = C19H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C9DCO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC9DCO2(i)=fC9DCO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H26O6(i)=fC19H26O6(i)+1; 

i=i+1;
Rnames{i} = 'C915O2 + APC10H15O4RB = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C915O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC915O2(i)=fC915O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C917O2 + APC10H15O4RB = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C917O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC917O2(i)=fC917O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C919O2 + APC10H15O4RB = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C919O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC919O2(i)=fC919O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C914O2 + APC10H15O4RB = C19H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C914O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC914O2(i)=fC914O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H26O7(i)=fC19H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C916O2 + APC10H15O4RB = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C916O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC916O2(i)=fC916O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C88CO3 + APC10H15O4RB = C19H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC88CO3(i)=fC88CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H26O7(i)=fC19H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C87CO3 + APC10H15O4RB = C19H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C87CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC87CO3(i)=fC87CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H26O8(i)=fC19H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C822CO3 + APC10H15O4RB = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C822CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC822CO3(i)=fC822CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'NLMKAO2 + APC10H15O4RB = C19H29O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLMKAO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNLMKAO2(i)=fNLMKAO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC19H29O5NO3(i)=fC19H29O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'C85O2 + APC10H15O4RB = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC85O2(i)=fC85O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O5(i)=fC18H28O5(i)+1; 

i=i+1;
Rnames{i} = 'C89O2 + APC10H15O4RB = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC89O2(i)=fC89O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O5(i)=fC18H28O5(i)+1; 

i=i+1;
Rnames{i} = 'C86O2 + APC10H15O4RB = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC86O2(i)=fC86O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C811O2 + APC10H15O4RB = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC811O2(i)=fC811O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C810O2 + APC10H15O4RB = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC810O2(i)=fC810O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C812O2 + APC10H15O4RB = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC812O2(i)=fC812O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C813O2 + APC10H15O4RB = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC813O2(i)=fC813O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O5(i)=fC18H28O5(i)+1; 

i=i+1;
Rnames{i} = 'C729CO3 + APC10H15O4RB = C18H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C729CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC729CO3(i)=fC729CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H26O6(i)=fC18H26O6(i)+1; 

i=i+1;
Rnames{i} = 'C816O2 + APC10H15O4RB = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C816O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC816O2(i)=fC816O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O5(i)=fC18H28O5(i)+1; 

i=i+1;
Rnames{i} = 'C817O2 + APC10H15O4RB = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C817O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC817O2(i)=fC817O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C826O2 + APC10H15O4RB = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C826O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC826O2(i)=fC826O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C822O2 + APC10H15O4RB = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C822O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC822O2(i)=fC822O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O5(i)=fC18H28O5(i)+1; 

i=i+1;
Rnames{i} = 'C818O2 + APC10H15O4RB = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C818O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC818O2(i)=fC818O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C823O2 + APC10H15O4RB = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C823O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC823O2(i)=fC823O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C819O2 + APC10H15O4RB = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C819O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC819O2(i)=fC819O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C727CO3 + APC10H15O4RB = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C727CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC727CO3(i)=fC727CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C731CO3 + APC10H15O4RB = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C731CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC731CO3(i)=fC731CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C824O2 + APC10H15O4RB = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C824O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC824O2(i)=fC824O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C820O2 + APC10H15O4RB = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C820O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC820O2(i)=fC820O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C825O2 + APC10H15O4RB = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C825O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC825O2(i)=fC825O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C732CO3 + APC10H15O4RB = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C732CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC732CO3(i)=fC732CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C8BCO2 + APC10H15O4RB = C18H28O4';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C8BCO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC8BCO2(i)=fC8BCO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H28O4(i)=fC18H28O4(i)+1; 

i=i+1;
Rnames{i} = 'C88O2 + APC10H15O4RB = C18H26O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C88O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC88O2(i)=fC88O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H26O6(i)=fC18H26O6(i)+1; 

i=i+1;
Rnames{i} = 'C718CO3 + APC10H15O4RB = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C718CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC718CO3(i)=fC718CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C87O2 + APC10H15O4RB = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C87O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fC87O2(i)=fC87O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C721CO3 + APC10H15O4RB = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'APC10H15O4RB'; 
fC721CO3(i)=fC721CO3(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'NC826O2 + APC10H15O4RB = C18H27O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NC826O2'; Gstr{i,2} = 'APC10H15O4RB'; 
fNC826O2(i)=fNC826O2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fC18H27O6NO3(i)=fC18H27O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'BUT2OLO2 + APC10H15O4RB = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'BUT2OLO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fBUT2OLO2(i)=fBUT2OLO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1; 

i=i+1;
Rnames{i} = 'CHEXO2 + APC10H15O4RB = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'CHEXO2'; Gstr{i,2} = 'APC10H15O4RB'; 
fCHEXO2(i)=fCHEXO2(i)-1; fAPC10H15O4RB(i)=fAPC10H15O4RB(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1; 

%% PINAL4H1RO2 rxns
i=i+1;
Rnames{i} = 'APINAO2 + PINAL4H1RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'APINBO2 + PINAL4H1RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'APINCO2 + PINAL4H1RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'C107O2 + PINAL4H1RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC107O2(i)=fC107O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C109O2 + PINAL4H1RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC109O2(i)=fC109O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C106O2 + PINAL4H1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC106O2(i)=fC106O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C920CO3 + PINAL4H1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC920CO3(i)=fC920CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C108O2 + PINAL4H1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC108O2(i)=fC108O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'PINALO2 + PINAL4H1RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fPINALO2(i)=fPINALO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C96CO3 + PINAL4H1RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC96CO3(i)=fC96CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C923CO3 + PINAL4H1RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C923CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC923CO3(i)=fC923CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMAO2 + PINAL4H1RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMAO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fLIMAO2(i)=fLIMAO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'LIMALBO2 + PINAL4H1RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALBO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fLIMALBO2(i)=fLIMALBO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMCO2 + PINAL4H1RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMCO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fLIMCO2(i)=fLIMCO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'LIMALO2 + PINAL4H1RO2 = C20H32O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fLIMALO2(i)=fLIMALO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H32O8(i)=fC20H32O8(i)+1; 

i=i+1;
Rnames{i} = 'LIMBO2 + PINAL4H1RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMBO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fLIMBO2(i)=fLIMBO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'LIMALAO2 + PINAL4H1RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALAO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fLIMALAO2(i)=fLIMALAO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'NAPINAO2 + PINAL4H1RO2 = C20H31O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINAO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNAPINAO2(i)=fNAPINAO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H31O5NO3(i)=fC20H31O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'NAPINBO2 + PINAL4H1RO2 = C20H31O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINBO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNAPINBO2(i)=fNAPINBO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H31O5NO3(i)=fC20H31O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'BPINAO2 + PINAL4H1RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'BPINAO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fBPINAO2(i)=fBPINAO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'BPINBO2 + PINAL4H1RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'BPINBO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fBPINBO2(i)=fBPINBO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'BPINCO2 + PINAL4H1RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'BPINCO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fBPINCO2(i)=fBPINCO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'C918CO3 + PINAL4H1RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC918CO3(i)=fC918CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'NC102O2 + PINAL4H1RO2 = C20H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC102O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNC102O2(i)=fNC102O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H29O7NO3(i)=fC20H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC101O2 + PINAL4H1RO2 = C20H29O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NC101O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNC101O2(i)=fNC101O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H29O6NO3(i)=fC20H29O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMO2 + PINAL4H1RO2 = C20H31O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLIMO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNLIMO2(i)=fNLIMO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H31O5NO3(i)=fC20H31O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMALO2 + PINAL4H1RO2 = C20H31O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLIMALO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNLIMALO2(i)=fNLIMALO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H31O7NO3(i)=fC20H31O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC91CO3 + PINAL4H1RO2 = C20H29O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NC91CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNC91CO3(i)=fNC91CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H29O6NO3(i)=fC20H29O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINAO2 + PINAL4H1RO2 = C20H31O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINAO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNBPINAO2(i)=fNBPINAO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H31O5NO3(i)=fC20H31O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINBO2 + PINAL4H1RO2 = C20H31O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINBO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNBPINBO2(i)=fNBPINBO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC20H31O5NO3(i)=fC20H31O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'C96O2 + PINAL4H1RO2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC96O2(i)=fC96O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C89CO3 + PINAL4H1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC89CO3(i)=fC89CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C920O2 + PINAL4H1RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC920O2(i)=fC920O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C97O2 + PINAL4H1RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC97O2(i)=fC97O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C85CO3 + PINAL4H1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC85CO3(i)=fC85CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C811CO3 + PINAL4H1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC811CO3(i)=fC811CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C921O2 + PINAL4H1RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC921O2(i)=fC921O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C98O2 + PINAL4H1RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC98O2(i)=fC98O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C922O2 + PINAL4H1RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC922O2(i)=fC922O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C923O2 + PINAL4H1RO2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C923O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC923O2(i)=fC923O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C924O2 + PINAL4H1RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C924O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC924O2(i)=fC924O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C816CO3 + PINAL4H1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C816CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC816CO3(i)=fC816CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NORLIMO2 + PINAL4H1RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NORLIMO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNORLIMO2(i)=fNORLIMO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LMKAO2 + PINAL4H1RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKAO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fLMKAO2(i)=fLMKAO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'LMKBO2 + PINAL4H1RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKBO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fLMKBO2(i)=fLMKBO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C926O2 + PINAL4H1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C926O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC926O2(i)=fC926O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C817CO3 + PINAL4H1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC817CO3(i)=fC817CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'LMLKAO2 + PINAL4H1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKAO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fLMLKAO2(i)=fLMLKAO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'LMLKBO2 + PINAL4H1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKBO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fLMLKBO2(i)=fLMLKBO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C823CO3 + PINAL4H1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC823CO3(i)=fC823CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C925O2 + PINAL4H1RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C925O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC925O2(i)=fC925O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'NOPINAO2 + PINAL4H1RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINAO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNOPINAO2(i)=fNOPINAO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'NOPINBO2 + PINAL4H1RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINBO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNOPINBO2(i)=fNOPINBO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'NOPINCO2 + PINAL4H1RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINCO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNOPINCO2(i)=fNOPINCO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'NOPINDO2 + PINAL4H1RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINDO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNOPINDO2(i)=fNOPINDO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C918O2 + PINAL4H1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC918O2(i)=fC918O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C9DCO2 + PINAL4H1RO2 = C19H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C9DCO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC9DCO2(i)=fC9DCO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H26O7(i)=fC19H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C915O2 + PINAL4H1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C915O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC915O2(i)=fC915O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C917O2 + PINAL4H1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C917O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC917O2(i)=fC917O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C919O2 + PINAL4H1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C919O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC919O2(i)=fC919O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C914O2 + PINAL4H1RO2 = C19H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C914O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC914O2(i)=fC914O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H26O8(i)=fC19H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C916O2 + PINAL4H1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C916O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC916O2(i)=fC916O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C88CO3 + PINAL4H1RO2 = C19H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC88CO3(i)=fC88CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H26O8(i)=fC19H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C87CO3 + PINAL4H1RO2 = C19H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C87CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC87CO3(i)=fC87CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H26O9(i)=fC19H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C822CO3 + PINAL4H1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C822CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC822CO3(i)=fC822CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NLMKAO2 + PINAL4H1RO2 = C19H29O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLMKAO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNLMKAO2(i)=fNLMKAO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC19H29O6NO3(i)=fC19H29O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'C85O2 + PINAL4H1RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC85O2(i)=fC85O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C89O2 + PINAL4H1RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC89O2(i)=fC89O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C86O2 + PINAL4H1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC86O2(i)=fC86O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C811O2 + PINAL4H1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC811O2(i)=fC811O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C810O2 + PINAL4H1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC810O2(i)=fC810O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C812O2 + PINAL4H1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC812O2(i)=fC812O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C813O2 + PINAL4H1RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC813O2(i)=fC813O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C729CO3 + PINAL4H1RO2 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C729CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC729CO3(i)=fC729CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C816O2 + PINAL4H1RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C816O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC816O2(i)=fC816O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C817O2 + PINAL4H1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC817O2(i)=fC817O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C826O2 + PINAL4H1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C826O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC826O2(i)=fC826O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C822O2 + PINAL4H1RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C822O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC822O2(i)=fC822O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C818O2 + PINAL4H1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C818O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC818O2(i)=fC818O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C823O2 + PINAL4H1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC823O2(i)=fC823O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C819O2 + PINAL4H1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C819O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC819O2(i)=fC819O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C727CO3 + PINAL4H1RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C727CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC727CO3(i)=fC727CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C731CO3 + PINAL4H1RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C731CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC731CO3(i)=fC731CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C824O2 + PINAL4H1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C824O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC824O2(i)=fC824O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C820O2 + PINAL4H1RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C820O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC820O2(i)=fC820O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C825O2 + PINAL4H1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C825O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC825O2(i)=fC825O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C732CO3 + PINAL4H1RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C732CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC732CO3(i)=fC732CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C8BCO2 + PINAL4H1RO2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C8BCO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC8BCO2(i)=fC8BCO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1; 

i=i+1;
Rnames{i} = 'C88O2 + PINAL4H1RO2 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC88O2(i)=fC88O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C718CO3 + PINAL4H1RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C718CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC718CO3(i)=fC718CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C87O2 + PINAL4H1RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C87O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC87O2(i)=fC87O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C721CO3 + PINAL4H1RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'PINAL4H1RO2'; 
fC721CO3(i)=fC721CO3(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'NC826O2 + PINAL4H1RO2 = C18H27O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC826O2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fNC826O2(i)=fNC826O2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fC18H27O7NO3(i)=fC18H27O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'BUT2OLO2 + PINAL4H1RO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'BUT2OLO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fBUT2OLO2(i)=fBUT2OLO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1; 

i=i+1;
Rnames{i} = 'CHEXO2 + PINAL4H1RO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'CHEXO2'; Gstr{i,2} = 'PINAL4H1RO2'; 
fCHEXO2(i)=fCHEXO2(i)-1; fPINAL4H1RO2(i)=fPINAL4H1RO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1; 

%% APC10H15O5RO2 rxns
i=i+1;
Rnames{i} = 'APINAO2 + APC10H15O5RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'APINBO2 + APC10H15O5RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'APINCO2 + APC10H15O5RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'C107O2 + APC10H15O5RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC107O2(i)=fC107O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C109O2 + APC10H15O5RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC109O2(i)=fC109O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C106O2 + APC10H15O5RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC106O2(i)=fC106O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C920CO3 + APC10H15O5RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC920CO3(i)=fC920CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C108O2 + APC10H15O5RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC108O2(i)=fC108O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'PINALO2 + APC10H15O5RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fPINALO2(i)=fPINALO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C96CO3 + APC10H15O5RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC96CO3(i)=fC96CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C923CO3 + APC10H15O5RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C923CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC923CO3(i)=fC923CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMAO2 + APC10H15O5RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMAO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fLIMAO2(i)=fLIMAO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'LIMALBO2 + APC10H15O5RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALBO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fLIMALBO2(i)=fLIMALBO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMCO2 + APC10H15O5RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMCO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fLIMCO2(i)=fLIMCO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'LIMALO2 + APC10H15O5RO2 = C20H32O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fLIMALO2(i)=fLIMALO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H32O8(i)=fC20H32O8(i)+1; 

i=i+1;
Rnames{i} = 'LIMBO2 + APC10H15O5RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'LIMBO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fLIMBO2(i)=fLIMBO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'LIMALAO2 + APC10H15O5RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALAO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fLIMALAO2(i)=fLIMALAO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'NAPINAO2 + APC10H15O5RO2 = C20H31O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINAO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNAPINAO2(i)=fNAPINAO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H31O5NO3(i)=fC20H31O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'NAPINBO2 + APC10H15O5RO2 = C20H31O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINBO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNAPINBO2(i)=fNAPINBO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H31O5NO3(i)=fC20H31O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'BPINAO2 + APC10H15O5RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'BPINAO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fBPINAO2(i)=fBPINAO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'BPINBO2 + APC10H15O5RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'BPINBO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fBPINBO2(i)=fBPINBO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'BPINCO2 + APC10H15O5RO2 = C20H32O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'BPINCO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fBPINCO2(i)=fBPINCO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H32O6(i)=fC20H32O6(i)+1; 

i=i+1;
Rnames{i} = 'C918CO3 + APC10H15O5RO2 = C20H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC918CO3(i)=fC918CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H30O7(i)=fC20H30O7(i)+1; 

i=i+1;
Rnames{i} = 'NC102O2 + APC10H15O5RO2 = C20H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC102O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNC102O2(i)=fNC102O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H29O7NO3(i)=fC20H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC101O2 + APC10H15O5RO2 = C20H29O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NC101O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNC101O2(i)=fNC101O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H29O6NO3(i)=fC20H29O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMO2 + APC10H15O5RO2 = C20H31O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLIMO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNLIMO2(i)=fNLIMO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H31O5NO3(i)=fC20H31O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMALO2 + APC10H15O5RO2 = C20H31O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLIMALO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNLIMALO2(i)=fNLIMALO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H31O7NO3(i)=fC20H31O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC91CO3 + APC10H15O5RO2 = C20H29O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NC91CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNC91CO3(i)=fNC91CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H29O6NO3(i)=fC20H29O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINAO2 + APC10H15O5RO2 = C20H31O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINAO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNBPINAO2(i)=fNBPINAO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H31O5NO3(i)=fC20H31O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINBO2 + APC10H15O5RO2 = C20H31O5NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINBO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNBPINBO2(i)=fNBPINBO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC20H31O5NO3(i)=fC20H31O5NO3(i)+1; 

i=i+1;
Rnames{i} = 'C96O2 + APC10H15O5RO2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC96O2(i)=fC96O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C89CO3 + APC10H15O5RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC89CO3(i)=fC89CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C920O2 + APC10H15O5RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC920O2(i)=fC920O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C97O2 + APC10H15O5RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC97O2(i)=fC97O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C85CO3 + APC10H15O5RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC85CO3(i)=fC85CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C811CO3 + APC10H15O5RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC811CO3(i)=fC811CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C921O2 + APC10H15O5RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC921O2(i)=fC921O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C98O2 + APC10H15O5RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC98O2(i)=fC98O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C922O2 + APC10H15O5RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC922O2(i)=fC922O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C923O2 + APC10H15O5RO2 = C19H30O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C923O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC923O2(i)=fC923O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H30O6(i)=fC19H30O6(i)+1; 

i=i+1;
Rnames{i} = 'C924O2 + APC10H15O5RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C924O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC924O2(i)=fC924O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C816CO3 + APC10H15O5RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C816CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC816CO3(i)=fC816CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NORLIMO2 + APC10H15O5RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NORLIMO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNORLIMO2(i)=fNORLIMO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LMKAO2 + APC10H15O5RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKAO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fLMKAO2(i)=fLMKAO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'LMKBO2 + APC10H15O5RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKBO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fLMKBO2(i)=fLMKBO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C926O2 + APC10H15O5RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C926O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC926O2(i)=fC926O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C817CO3 + APC10H15O5RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC817CO3(i)=fC817CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'LMLKAO2 + APC10H15O5RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKAO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fLMLKAO2(i)=fLMLKAO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'LMLKBO2 + APC10H15O5RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKBO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fLMLKBO2(i)=fLMLKBO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C823CO3 + APC10H15O5RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC823CO3(i)=fC823CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C925O2 + APC10H15O5RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C925O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC925O2(i)=fC925O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'NOPINAO2 + APC10H15O5RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINAO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNOPINAO2(i)=fNOPINAO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'NOPINBO2 + APC10H15O5RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINBO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNOPINBO2(i)=fNOPINBO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'NOPINCO2 + APC10H15O5RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINCO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNOPINCO2(i)=fNOPINCO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'NOPINDO2 + APC10H15O5RO2 = C19H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'NOPINDO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNOPINDO2(i)=fNOPINDO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O6(i)=fC19H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C918O2 + APC10H15O5RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC918O2(i)=fC918O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C9DCO2 + APC10H15O5RO2 = C19H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C9DCO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC9DCO2(i)=fC9DCO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H26O7(i)=fC19H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C915O2 + APC10H15O5RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C915O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC915O2(i)=fC915O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C917O2 + APC10H15O5RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C917O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC917O2(i)=fC917O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C919O2 + APC10H15O5RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C919O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC919O2(i)=fC919O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C914O2 + APC10H15O5RO2 = C19H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C914O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC914O2(i)=fC914O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H26O8(i)=fC19H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C916O2 + APC10H15O5RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C916O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC916O2(i)=fC916O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C88CO3 + APC10H15O5RO2 = C19H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC88CO3(i)=fC88CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H26O8(i)=fC19H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C87CO3 + APC10H15O5RO2 = C19H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C87CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC87CO3(i)=fC87CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H26O9(i)=fC19H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C822CO3 + APC10H15O5RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C822CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC822CO3(i)=fC822CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NLMKAO2 + APC10H15O5RO2 = C19H29O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLMKAO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNLMKAO2(i)=fNLMKAO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC19H29O6NO3(i)=fC19H29O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'C85O2 + APC10H15O5RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC85O2(i)=fC85O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C89O2 + APC10H15O5RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC89O2(i)=fC89O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C86O2 + APC10H15O5RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC86O2(i)=fC86O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C811O2 + APC10H15O5RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC811O2(i)=fC811O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C810O2 + APC10H15O5RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC810O2(i)=fC810O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C812O2 + APC10H15O5RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC812O2(i)=fC812O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H26O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C813O2 + APC10H15O5RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC813O2(i)=fC813O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C729CO3 + APC10H15O5RO2 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C729CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC729CO3(i)=fC729CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C816O2 + APC10H15O5RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C816O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC816O2(i)=fC816O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C817O2 + APC10H15O5RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC817O2(i)=fC817O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C826O2 + APC10H15O5RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C826O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC826O2(i)=fC826O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C822O2 + APC10H15O5RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C822O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC822O2(i)=fC822O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C818O2 + APC10H15O5RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C818O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC818O2(i)=fC818O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C823O2 + APC10H15O5RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC823O2(i)=fC823O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C819O2 + APC10H15O5RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C819O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC819O2(i)=fC819O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C727CO3 + APC10H15O5RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C727CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC727CO3(i)=fC727CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C731CO3 + APC10H15O5RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C731CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC731CO3(i)=fC731CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C824O2 + APC10H15O5RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C824O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC824O2(i)=fC824O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C820O2 + APC10H15O5RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C820O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC820O2(i)=fC820O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C825O2 + APC10H15O5RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C825O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC825O2(i)=fC825O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C732CO3 + APC10H15O5RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C732CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC732CO3(i)=fC732CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C8BCO2 + APC10H15O5RO2 = C18H28O5';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C8BCO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC8BCO2(i)=fC8BCO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H28O5(i)=fC18H28O5(i)+1; 

i=i+1;
Rnames{i} = 'C88O2 + APC10H15O5RO2 = C18H26O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC88O2(i)=fC88O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H26O7(i)=fC18H26O7(i)+1; 

i=i+1;
Rnames{i} = 'C718CO3 + APC10H15O5RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C718CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC718CO3(i)=fC718CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C87O2 + APC10H15O5RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C87O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC87O2(i)=fC87O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C721CO3 + APC10H15O5RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'APC10H15O5RO2'; 
fC721CO3(i)=fC721CO3(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'NC826O2 + APC10H15O5RO2 = C18H27O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC826O2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fNC826O2(i)=fNC826O2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fC18H27O7NO3(i)=fC18H27O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'BUT2OLO2 + APC10H15O5RO2 = BUT2OLdimer';
k(:,i) = kROORa.*fROOR1;
Gstr{i,1} = 'BUT2OLO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fBUT2OLO2(i)=fBUT2OLO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1; 

i=i+1;
Rnames{i} = 'CHEXO2 + APC10H15O5RO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'CHEXO2'; Gstr{i,2} = 'APC10H15O5RO2'; 
fCHEXO2(i)=fCHEXO2(i)-1; fAPC10H15O5RO2(i)=fAPC10H15O5RO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1; 

%% PINALPA4RO2 rxns
i=i+1;
Rnames{i} = 'APINAO2 + PINALPA4RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'APINBO2 + PINALPA4RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'APINCO2 + PINALPA4RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'C107O2 + PINALPA4RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC107O2(i)=fC107O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C109O2 + PINALPA4RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC109O2(i)=fC109O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C106O2 + PINALPA4RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC106O2(i)=fC106O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C920CO3 + PINALPA4RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC920CO3(i)=fC920CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C108O2 + PINALPA4RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC108O2(i)=fC108O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'PINALO2 + PINALPA4RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fPINALO2(i)=fPINALO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C96CO3 + PINALPA4RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC96CO3(i)=fC96CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C923CO3 + PINALPA4RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C923CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC923CO3(i)=fC923CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LIMAO2 + PINALPA4RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMAO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fLIMAO2(i)=fLIMAO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMALBO2 + PINALPA4RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALBO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fLIMALBO2(i)=fLIMALBO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LIMCO2 + PINALPA4RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMCO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fLIMCO2(i)=fLIMCO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMALO2 + PINALPA4RO2 = C20H32O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fLIMALO2(i)=fLIMALO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H32O9(i)=fC20H32O9(i)+1; 

i=i+1;
Rnames{i} = 'LIMBO2 + PINALPA4RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMBO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fLIMBO2(i)=fLIMBO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMALAO2 + PINALPA4RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALAO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fLIMALAO2(i)=fLIMALAO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'NAPINAO2 + PINALPA4RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINAO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNAPINAO2(i)=fNAPINAO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NAPINBO2 + PINALPA4RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINBO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNAPINBO2(i)=fNAPINBO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'BPINAO2 + PINALPA4RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINAO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fBPINAO2(i)=fBPINAO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'BPINBO2 + PINALPA4RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINBO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fBPINBO2(i)=fBPINBO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'BPINCO2 + PINALPA4RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINCO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fBPINCO2(i)=fBPINCO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'C918CO3 + PINALPA4RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC918CO3(i)=fC918CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'NC102O2 + PINALPA4RO2 = C20H29O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC102O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNC102O2(i)=fNC102O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H29O8NO3(i)=fC20H29O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC101O2 + PINALPA4RO2 = C20H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC101O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNC101O2(i)=fNC101O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H29O7NO3(i)=fC20H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMO2 + PINALPA4RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLIMO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNLIMO2(i)=fNLIMO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMALO2 + PINALPA4RO2 = C20H31O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLIMALO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNLIMALO2(i)=fNLIMALO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H31O8NO3(i)=fC20H31O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC91CO3 + PINALPA4RO2 = C20H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC91CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fNC91CO3(i)=fNC91CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H29O7NO3(i)=fC20H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINAO2 + PINALPA4RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINAO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNBPINAO2(i)=fNBPINAO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINBO2 + PINALPA4RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINBO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNBPINBO2(i)=fNBPINBO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'C96O2 + PINALPA4RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC96O2(i)=fC96O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C89CO3 + PINALPA4RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC89CO3(i)=fC89CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C920O2 + PINALPA4RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC920O2(i)=fC920O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C97O2 + PINALPA4RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC97O2(i)=fC97O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C85CO3 + PINALPA4RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC85CO3(i)=fC85CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C811CO3 + PINALPA4RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC811CO3(i)=fC811CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C921O2 + PINALPA4RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC921O2(i)=fC921O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C98O2 + PINALPA4RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC98O2(i)=fC98O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C922O2 + PINALPA4RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC922O2(i)=fC922O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C923O2 + PINALPA4RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C923O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC923O2(i)=fC923O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C924O2 + PINALPA4RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C924O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC924O2(i)=fC924O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C816CO3 + PINALPA4RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C816CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC816CO3(i)=fC816CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'NORLIMO2 + PINALPA4RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NORLIMO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNORLIMO2(i)=fNORLIMO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'LMKAO2 + PINALPA4RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKAO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fLMKAO2(i)=fLMKAO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LMKBO2 + PINALPA4RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKBO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fLMKBO2(i)=fLMKBO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C926O2 + PINALPA4RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C926O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC926O2(i)=fC926O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C817CO3 + PINALPA4RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC817CO3(i)=fC817CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'LMLKAO2 + PINALPA4RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKAO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fLMLKAO2(i)=fLMLKAO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'LMLKBO2 + PINALPA4RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKBO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fLMLKBO2(i)=fLMLKBO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C823CO3 + PINALPA4RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC823CO3(i)=fC823CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C925O2 + PINALPA4RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C925O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC925O2(i)=fC925O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'NOPINAO2 + PINALPA4RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINAO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNOPINAO2(i)=fNOPINAO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NOPINBO2 + PINALPA4RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINBO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNOPINBO2(i)=fNOPINBO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NOPINCO2 + PINALPA4RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINCO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNOPINCO2(i)=fNOPINCO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NOPINDO2 + PINALPA4RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINDO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNOPINDO2(i)=fNOPINDO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C918O2 + PINALPA4RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC918O2(i)=fC918O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C9DCO2 + PINALPA4RO2 = C19H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C9DCO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC9DCO2(i)=fC9DCO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H26O8(i)=fC19H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C915O2 + PINALPA4RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C915O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC915O2(i)=fC915O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C917O2 + PINALPA4RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C917O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC917O2(i)=fC917O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C919O2 + PINALPA4RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C919O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC919O2(i)=fC919O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C914O2 + PINALPA4RO2 = C19H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C914O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC914O2(i)=fC914O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H26O9(i)=fC19H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C916O2 + PINALPA4RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C916O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC916O2(i)=fC916O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C88CO3 + PINALPA4RO2 = C19H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC88CO3(i)=fC88CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H26O9(i)=fC19H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C87CO3 + PINALPA4RO2 = C19H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C87CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC87CO3(i)=fC87CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H26O10(i)=fC19H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C822CO3 + PINALPA4RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C822CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC822CO3(i)=fC822CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'NLMKAO2 + PINALPA4RO2 = C19H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLMKAO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNLMKAO2(i)=fNLMKAO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC19H29O7NO3(i)=fC19H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'C85O2 + PINALPA4RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC85O2(i)=fC85O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C89O2 + PINALPA4RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC89O2(i)=fC89O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C86O2 + PINALPA4RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC86O2(i)=fC86O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C811O2 + PINALPA4RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC811O2(i)=fC811O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C810O2 + PINALPA4RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC810O2(i)=fC810O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C812O2 + PINALPA4RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC812O2(i)=fC812O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C813O2 + PINALPA4RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC813O2(i)=fC813O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C729CO3 + PINALPA4RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C729CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC729CO3(i)=fC729CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C816O2 + PINALPA4RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C816O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC816O2(i)=fC816O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C817O2 + PINALPA4RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC817O2(i)=fC817O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C826O2 + PINALPA4RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C826O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC826O2(i)=fC826O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C822O2 + PINALPA4RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C822O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC822O2(i)=fC822O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C818O2 + PINALPA4RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C818O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC818O2(i)=fC818O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C823O2 + PINALPA4RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC823O2(i)=fC823O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C819O2 + PINALPA4RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C819O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC819O2(i)=fC819O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C727CO3 + PINALPA4RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C727CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC727CO3(i)=fC727CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C731CO3 + PINALPA4RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C731CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC731CO3(i)=fC731CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C824O2 + PINALPA4RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C824O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC824O2(i)=fC824O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C820O2 + PINALPA4RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C820O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC820O2(i)=fC820O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C825O2 + PINALPA4RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C825O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC825O2(i)=fC825O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C732CO3 + PINALPA4RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C732CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC732CO3(i)=fC732CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C8BCO2 + PINALPA4RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C8BCO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC8BCO2(i)=fC8BCO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C88O2 + PINALPA4RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC88O2(i)=fC88O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C718CO3 + PINALPA4RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C718CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC718CO3(i)=fC718CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C87O2 + PINALPA4RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C87O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fC87O2(i)=fC87O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C721CO3 + PINALPA4RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'PINALPA4RO2'; 
fC721CO3(i)=fC721CO3(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'NC826O2 + PINALPA4RO2 = C18H27O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC826O2'; Gstr{i,2} = 'PINALPA4RO2'; 
fNC826O2(i)=fNC826O2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fC18H27O8NO3(i)=fC18H27O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'BUT2OLO2 + PINALPA4RO2 = BUT2OLdimer';
k(:,i) = kROORb.*fROOR1;
Gstr{i,1} = 'BUT2OLO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fBUT2OLO2(i)=fBUT2OLO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1; 

i=i+1;
Rnames{i} = 'CHEXO2 + PINALPA4RO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'CHEXO2'; Gstr{i,2} = 'PINALPA4RO2'; 
fCHEXO2(i)=fCHEXO2(i)-1; fPINALPA4RO2(i)=fPINALPA4RO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1; 

%% APC10H15O6RO2 rxns
i=i+1;
Rnames{i} = 'APINAO2 + APC10H15O6RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'APINBO2 + APC10H15O6RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'APINCO2 + APC10H15O6RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'C107O2 + APC10H15O6RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC107O2(i)=fC107O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C109O2 + APC10H15O6RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC109O2(i)=fC109O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C106O2 + APC10H15O6RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC106O2(i)=fC106O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C920CO3 + APC10H15O6RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC920CO3(i)=fC920CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C108O2 + APC10H15O6RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC108O2(i)=fC108O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'PINALO2 + APC10H15O6RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fPINALO2(i)=fPINALO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C96CO3 + APC10H15O6RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC96CO3(i)=fC96CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C923CO3 + APC10H15O6RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C923CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC923CO3(i)=fC923CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LIMAO2 + APC10H15O6RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMAO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fLIMAO2(i)=fLIMAO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMALBO2 + APC10H15O6RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALBO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fLIMALBO2(i)=fLIMALBO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LIMCO2 + APC10H15O6RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMCO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fLIMCO2(i)=fLIMCO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMALO2 + APC10H15O6RO2 = C20H32O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fLIMALO2(i)=fLIMALO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H32O9(i)=fC20H32O9(i)+1; 

i=i+1;
Rnames{i} = 'LIMBO2 + APC10H15O6RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMBO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fLIMBO2(i)=fLIMBO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMALAO2 + APC10H15O6RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALAO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fLIMALAO2(i)=fLIMALAO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'NAPINAO2 + APC10H15O6RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINAO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNAPINAO2(i)=fNAPINAO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NAPINBO2 + APC10H15O6RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINBO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNAPINBO2(i)=fNAPINBO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'BPINAO2 + APC10H15O6RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINAO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fBPINAO2(i)=fBPINAO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'BPINBO2 + APC10H15O6RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINBO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fBPINBO2(i)=fBPINBO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'BPINCO2 + APC10H15O6RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINCO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fBPINCO2(i)=fBPINCO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'C918CO3 + APC10H15O6RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC918CO3(i)=fC918CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'NC102O2 + APC10H15O6RO2 = C20H29O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC102O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNC102O2(i)=fNC102O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H29O8NO3(i)=fC20H29O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC101O2 + APC10H15O6RO2 = C20H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC101O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNC101O2(i)=fNC101O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H29O7NO3(i)=fC20H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMO2 + APC10H15O6RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLIMO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNLIMO2(i)=fNLIMO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMALO2 + APC10H15O6RO2 = C20H31O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLIMALO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNLIMALO2(i)=fNLIMALO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H31O8NO3(i)=fC20H31O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC91CO3 + APC10H15O6RO2 = C20H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC91CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNC91CO3(i)=fNC91CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H29O7NO3(i)=fC20H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINAO2 + APC10H15O6RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINAO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNBPINAO2(i)=fNBPINAO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINBO2 + APC10H15O6RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINBO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNBPINBO2(i)=fNBPINBO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'C96O2 + APC10H15O6RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC96O2(i)=fC96O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C89CO3 + APC10H15O6RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC89CO3(i)=fC89CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C920O2 + APC10H15O6RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC920O2(i)=fC920O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C97O2 + APC10H15O6RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC97O2(i)=fC97O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C85CO3 + APC10H15O6RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC85CO3(i)=fC85CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C811CO3 + APC10H15O6RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC811CO3(i)=fC811CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C921O2 + APC10H15O6RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC921O2(i)=fC921O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C98O2 + APC10H15O6RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC98O2(i)=fC98O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C922O2 + APC10H15O6RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC922O2(i)=fC922O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C923O2 + APC10H15O6RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C923O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC923O2(i)=fC923O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C924O2 + APC10H15O6RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C924O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC924O2(i)=fC924O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C816CO3 + APC10H15O6RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C816CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC816CO3(i)=fC816CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'NORLIMO2 + APC10H15O6RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NORLIMO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNORLIMO2(i)=fNORLIMO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'LMKAO2 + APC10H15O6RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKAO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fLMKAO2(i)=fLMKAO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LMKBO2 + APC10H15O6RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKBO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fLMKBO2(i)=fLMKBO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C926O2 + APC10H15O6RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C926O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC926O2(i)=fC926O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C817CO3 + APC10H15O6RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC817CO3(i)=fC817CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'LMLKAO2 + APC10H15O6RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKAO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fLMLKAO2(i)=fLMLKAO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'LMLKBO2 + APC10H15O6RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKBO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fLMLKBO2(i)=fLMLKBO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C823CO3 + APC10H15O6RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC823CO3(i)=fC823CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C925O2 + APC10H15O6RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C925O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC925O2(i)=fC925O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'NOPINAO2 + APC10H15O6RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINAO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNOPINAO2(i)=fNOPINAO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NOPINBO2 + APC10H15O6RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINBO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNOPINBO2(i)=fNOPINBO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NOPINCO2 + APC10H15O6RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINCO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNOPINCO2(i)=fNOPINCO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NOPINDO2 + APC10H15O6RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINDO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNOPINDO2(i)=fNOPINDO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C918O2 + APC10H15O6RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC918O2(i)=fC918O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C9DCO2 + APC10H15O6RO2 = C19H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C9DCO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC9DCO2(i)=fC9DCO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H26O8(i)=fC19H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C915O2 + APC10H15O6RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C915O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC915O2(i)=fC915O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C917O2 + APC10H15O6RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C917O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC917O2(i)=fC917O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C919O2 + APC10H15O6RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C919O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC919O2(i)=fC919O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C914O2 + APC10H15O6RO2 = C19H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C914O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC914O2(i)=fC914O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H26O9(i)=fC19H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C916O2 + APC10H15O6RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C916O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC916O2(i)=fC916O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C88CO3 + APC10H15O6RO2 = C19H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC88CO3(i)=fC88CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H26O9(i)=fC19H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C87CO3 + APC10H15O6RO2 = C19H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C87CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC87CO3(i)=fC87CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H26O10(i)=fC19H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C822CO3 + APC10H15O6RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C822CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC822CO3(i)=fC822CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'NLMKAO2 + APC10H15O6RO2 = C19H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLMKAO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNLMKAO2(i)=fNLMKAO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC19H29O7NO3(i)=fC19H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'C85O2 + APC10H15O6RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC85O2(i)=fC85O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C89O2 + APC10H15O6RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC89O2(i)=fC89O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C86O2 + APC10H15O6RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC86O2(i)=fC86O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C811O2 + APC10H15O6RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC811O2(i)=fC811O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C810O2 + APC10H15O6RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC810O2(i)=fC810O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C812O2 + APC10H15O6RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC812O2(i)=fC812O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C813O2 + APC10H15O6RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC813O2(i)=fC813O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C729CO3 + APC10H15O6RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C729CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC729CO3(i)=fC729CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C816O2 + APC10H15O6RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C816O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC816O2(i)=fC816O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C817O2 + APC10H15O6RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC817O2(i)=fC817O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C826O2 + APC10H15O6RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C826O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC826O2(i)=fC826O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C822O2 + APC10H15O6RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C822O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC822O2(i)=fC822O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C818O2 + APC10H15O6RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C818O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC818O2(i)=fC818O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C823O2 + APC10H15O6RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC823O2(i)=fC823O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C819O2 + APC10H15O6RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C819O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC819O2(i)=fC819O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C727CO3 + APC10H15O6RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C727CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC727CO3(i)=fC727CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C731CO3 + APC10H15O6RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C731CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC731CO3(i)=fC731CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C824O2 + APC10H15O6RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C824O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC824O2(i)=fC824O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C820O2 + APC10H15O6RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C820O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC820O2(i)=fC820O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C825O2 + APC10H15O6RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C825O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC825O2(i)=fC825O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C732CO3 + APC10H15O6RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C732CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC732CO3(i)=fC732CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C8BCO2 + APC10H15O6RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C8BCO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC8BCO2(i)=fC8BCO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C88O2 + APC10H15O6RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC88O2(i)=fC88O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C718CO3 + APC10H15O6RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C718CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC718CO3(i)=fC718CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C87O2 + APC10H15O6RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C87O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC87O2(i)=fC87O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C721CO3 + APC10H15O6RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'APC10H15O6RO2'; 
fC721CO3(i)=fC721CO3(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'NC826O2 + APC10H15O6RO2 = C18H27O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC826O2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fNC826O2(i)=fNC826O2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fC18H27O8NO3(i)=fC18H27O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'BUT2OLO2 + APC10H15O6RO2 = BUT2OLdimer';
k(:,i) = kROORb.*fROOR1;
Gstr{i,1} = 'BUT2OLO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fBUT2OLO2(i)=fBUT2OLO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1; 

i=i+1;
Rnames{i} = 'CHEXO2 + APC10H15O6RO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'CHEXO2'; Gstr{i,2} = 'APC10H15O6RO2'; 
fCHEXO2(i)=fCHEXO2(i)-1; fAPC10H15O6RO2(i)=fAPC10H15O6RO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1; 

%% PINAL10P1RO2 rxns
i=i+1;
Rnames{i} = 'APINAO2 + PINAL10P1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'APINBO2 + PINAL10P1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'APINCO2 + PINAL10P1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'C107O2 + PINAL10P1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC107O2(i)=fC107O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C109O2 + PINAL10P1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC109O2(i)=fC109O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C106O2 + PINAL10P1RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC106O2(i)=fC106O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C920CO3 + PINAL10P1RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC920CO3(i)=fC920CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C108O2 + PINAL10P1RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC108O2(i)=fC108O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'PINALO2 + PINAL10P1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fPINALO2(i)=fPINALO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C96CO3 + PINAL10P1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC96CO3(i)=fC96CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C923CO3 + PINAL10P1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C923CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC923CO3(i)=fC923CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LIMAO2 + PINAL10P1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMAO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fLIMAO2(i)=fLIMAO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMALBO2 + PINAL10P1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALBO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fLIMALBO2(i)=fLIMALBO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LIMCO2 + PINAL10P1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMCO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fLIMCO2(i)=fLIMCO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMALO2 + PINAL10P1RO2 = C20H32O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fLIMALO2(i)=fLIMALO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H32O9(i)=fC20H32O9(i)+1; 

i=i+1;
Rnames{i} = 'LIMBO2 + PINAL10P1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMBO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fLIMBO2(i)=fLIMBO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMALAO2 + PINAL10P1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALAO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fLIMALAO2(i)=fLIMALAO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'NAPINAO2 + PINAL10P1RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINAO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNAPINAO2(i)=fNAPINAO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NAPINBO2 + PINAL10P1RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINBO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNAPINBO2(i)=fNAPINBO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'BPINAO2 + PINAL10P1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINAO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fBPINAO2(i)=fBPINAO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'BPINBO2 + PINAL10P1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINBO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fBPINBO2(i)=fBPINBO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'BPINCO2 + PINAL10P1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINCO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fBPINCO2(i)=fBPINCO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'C918CO3 + PINAL10P1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC918CO3(i)=fC918CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'NC102O2 + PINAL10P1RO2 = C20H29O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC102O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNC102O2(i)=fNC102O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H29O8NO3(i)=fC20H29O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC101O2 + PINAL10P1RO2 = C20H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC101O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNC101O2(i)=fNC101O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H29O7NO3(i)=fC20H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMO2 + PINAL10P1RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLIMO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNLIMO2(i)=fNLIMO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMALO2 + PINAL10P1RO2 = C20H31O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLIMALO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNLIMALO2(i)=fNLIMALO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H31O8NO3(i)=fC20H31O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC91CO3 + PINAL10P1RO2 = C20H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC91CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNC91CO3(i)=fNC91CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H29O7NO3(i)=fC20H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINAO2 + PINAL10P1RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINAO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNBPINAO2(i)=fNBPINAO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINBO2 + PINAL10P1RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINBO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNBPINBO2(i)=fNBPINBO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'C96O2 + PINAL10P1RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC96O2(i)=fC96O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C89CO3 + PINAL10P1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC89CO3(i)=fC89CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C920O2 + PINAL10P1RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC920O2(i)=fC920O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C97O2 + PINAL10P1RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC97O2(i)=fC97O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C85CO3 + PINAL10P1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC85CO3(i)=fC85CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C811CO3 + PINAL10P1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC811CO3(i)=fC811CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C921O2 + PINAL10P1RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC921O2(i)=fC921O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C98O2 + PINAL10P1RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC98O2(i)=fC98O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C922O2 + PINAL10P1RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC922O2(i)=fC922O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C923O2 + PINAL10P1RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C923O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC923O2(i)=fC923O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C924O2 + PINAL10P1RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C924O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC924O2(i)=fC924O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C816CO3 + PINAL10P1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C816CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC816CO3(i)=fC816CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'NORLIMO2 + PINAL10P1RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NORLIMO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNORLIMO2(i)=fNORLIMO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'LMKAO2 + PINAL10P1RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKAO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fLMKAO2(i)=fLMKAO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LMKBO2 + PINAL10P1RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKBO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fLMKBO2(i)=fLMKBO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C926O2 + PINAL10P1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C926O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC926O2(i)=fC926O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C817CO3 + PINAL10P1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC817CO3(i)=fC817CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'LMLKAO2 + PINAL10P1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKAO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fLMLKAO2(i)=fLMLKAO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'LMLKBO2 + PINAL10P1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKBO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fLMLKBO2(i)=fLMLKBO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C823CO3 + PINAL10P1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC823CO3(i)=fC823CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C925O2 + PINAL10P1RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C925O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC925O2(i)=fC925O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'NOPINAO2 + PINAL10P1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINAO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNOPINAO2(i)=fNOPINAO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NOPINBO2 + PINAL10P1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINBO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNOPINBO2(i)=fNOPINBO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NOPINCO2 + PINAL10P1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINCO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNOPINCO2(i)=fNOPINCO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NOPINDO2 + PINAL10P1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINDO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNOPINDO2(i)=fNOPINDO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C918O2 + PINAL10P1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC918O2(i)=fC918O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C9DCO2 + PINAL10P1RO2 = C19H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C9DCO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC9DCO2(i)=fC9DCO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H26O8(i)=fC19H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C915O2 + PINAL10P1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C915O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC915O2(i)=fC915O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C917O2 + PINAL10P1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C917O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC917O2(i)=fC917O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C919O2 + PINAL10P1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C919O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC919O2(i)=fC919O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C914O2 + PINAL10P1RO2 = C19H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C914O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC914O2(i)=fC914O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H26O9(i)=fC19H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C916O2 + PINAL10P1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C916O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC916O2(i)=fC916O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C88CO3 + PINAL10P1RO2 = C19H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC88CO3(i)=fC88CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H26O9(i)=fC19H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C87CO3 + PINAL10P1RO2 = C19H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C87CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC87CO3(i)=fC87CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H26O10(i)=fC19H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C822CO3 + PINAL10P1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C822CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC822CO3(i)=fC822CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'NLMKAO2 + PINAL10P1RO2 = C19H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLMKAO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNLMKAO2(i)=fNLMKAO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC19H29O7NO3(i)=fC19H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'C85O2 + PINAL10P1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC85O2(i)=fC85O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C89O2 + PINAL10P1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC89O2(i)=fC89O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C86O2 + PINAL10P1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC86O2(i)=fC86O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C811O2 + PINAL10P1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC811O2(i)=fC811O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C810O2 + PINAL10P1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC810O2(i)=fC810O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C812O2 + PINAL10P1RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC812O2(i)=fC812O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C813O2 + PINAL10P1RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC813O2(i)=fC813O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C729CO3 + PINAL10P1RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C729CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC729CO3(i)=fC729CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C816O2 + PINAL10P1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C816O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC816O2(i)=fC816O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C817O2 + PINAL10P1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC817O2(i)=fC817O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C826O2 + PINAL10P1RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C826O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC826O2(i)=fC826O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C822O2 + PINAL10P1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C822O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC822O2(i)=fC822O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C818O2 + PINAL10P1RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C818O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC818O2(i)=fC818O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C823O2 + PINAL10P1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC823O2(i)=fC823O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C819O2 + PINAL10P1RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C819O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC819O2(i)=fC819O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C727CO3 + PINAL10P1RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C727CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC727CO3(i)=fC727CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C731CO3 + PINAL10P1RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C731CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC731CO3(i)=fC731CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C824O2 + PINAL10P1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C824O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC824O2(i)=fC824O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C820O2 + PINAL10P1RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C820O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC820O2(i)=fC820O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C825O2 + PINAL10P1RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C825O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC825O2(i)=fC825O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C732CO3 + PINAL10P1RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C732CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC732CO3(i)=fC732CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C8BCO2 + PINAL10P1RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C8BCO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC8BCO2(i)=fC8BCO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C88O2 + PINAL10P1RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC88O2(i)=fC88O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C718CO3 + PINAL10P1RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C718CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC718CO3(i)=fC718CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C87O2 + PINAL10P1RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C87O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC87O2(i)=fC87O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C721CO3 + PINAL10P1RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'PINAL10P1RO2'; 
fC721CO3(i)=fC721CO3(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'NC826O2 + PINAL10P1RO2 = C18H27O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC826O2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fNC826O2(i)=fNC826O2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fC18H27O8NO3(i)=fC18H27O8NO3(i)+1;

i=i+1;
Rnames{i} = 'BUT2OLO2 + PINAL10P1RO2 = BUT2OLdimer';
k(:,i) = kROORb.*fROOR1;
Gstr{i,1} = 'BUT2OLO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fBUT2OLO2(i)=fBUT2OLO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1; 

i=i+1;
Rnames{i} = 'CHEXO2 + PINAL10P1RO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'CHEXO2'; Gstr{i,2} = 'PINAL10P1RO2'; 
fCHEXO2(i)=fCHEXO2(i)-1; fPINAL10P1RO2(i)=fPINAL10P1RO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1; 

%% PINALPA1RO2 rxns
i=i+1;
Rnames{i} = 'APINAO2 + PINALPA1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'APINBO2 + PINALPA1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'APINCO2 + PINALPA1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'C107O2 + PINALPA1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC107O2(i)=fC107O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C109O2 + PINALPA1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC109O2(i)=fC109O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C106O2 + PINALPA1RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC106O2(i)=fC106O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C920CO3 + PINALPA1RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC920CO3(i)=fC920CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C108O2 + PINALPA1RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC108O2(i)=fC108O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'PINALO2 + PINALPA1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fPINALO2(i)=fPINALO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C96CO3 + PINALPA1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC96CO3(i)=fC96CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C923CO3 + PINALPA1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C923CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC923CO3(i)=fC923CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LIMAO2 + PINALPA1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMAO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fLIMAO2(i)=fLIMAO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMALBO2 + PINALPA1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALBO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fLIMALBO2(i)=fLIMALBO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LIMCO2 + PINALPA1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMCO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fLIMCO2(i)=fLIMCO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMALO2 + PINALPA1RO2 = C20H32O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fLIMALO2(i)=fLIMALO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H32O9(i)=fC20H32O9(i)+1; 

i=i+1;
Rnames{i} = 'LIMBO2 + PINALPA1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMBO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fLIMBO2(i)=fLIMBO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'LIMALAO2 + PINALPA1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALAO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fLIMALAO2(i)=fLIMALAO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'NAPINAO2 + PINALPA1RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINAO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNAPINAO2(i)=fNAPINAO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NAPINBO2 + PINALPA1RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NAPINBO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNAPINBO2(i)=fNAPINBO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'BPINAO2 + PINALPA1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINAO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fBPINAO2(i)=fBPINAO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'BPINBO2 + PINALPA1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINBO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fBPINBO2(i)=fBPINBO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'BPINCO2 + PINALPA1RO2 = C20H32O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINCO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fBPINCO2(i)=fBPINCO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H32O7(i)=fC20H32O7(i)+1; 

i=i+1;
Rnames{i} = 'C918CO3 + PINALPA1RO2 = C20H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC918CO3(i)=fC918CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H30O8(i)=fC20H30O8(i)+1; 

i=i+1;
Rnames{i} = 'NC102O2 + PINALPA1RO2 = C20H29O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC102O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNC102O2(i)=fNC102O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H29O8NO3(i)=fC20H29O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC101O2 + PINALPA1RO2 = C20H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC101O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNC101O2(i)=fNC101O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H29O7NO3(i)=fC20H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMO2 + PINALPA1RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NLIMO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNLIMO2(i)=fNLIMO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMALO2 + PINALPA1RO2 = C20H31O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLIMALO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNLIMALO2(i)=fNLIMALO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H31O8NO3(i)=fC20H31O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC91CO3 + PINALPA1RO2 = C20H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC91CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fNC91CO3(i)=fNC91CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H29O7NO3(i)=fC20H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINAO2 + PINALPA1RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINAO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNBPINAO2(i)=fNBPINAO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINBO2 + PINALPA1RO2 = C20H31O6NO3';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NBPINBO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNBPINBO2(i)=fNBPINBO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC20H31O6NO3(i)=fC20H31O6NO3(i)+1; 

i=i+1;
Rnames{i} = 'C96O2 + PINALPA1RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC96O2(i)=fC96O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C89CO3 + PINALPA1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC89CO3(i)=fC89CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C920O2 + PINALPA1RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC920O2(i)=fC920O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C97O2 + PINALPA1RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC97O2(i)=fC97O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C85CO3 + PINALPA1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC85CO3(i)=fC85CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C811CO3 + PINALPA1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC811CO3(i)=fC811CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C921O2 + PINALPA1RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC921O2(i)=fC921O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C98O2 + PINALPA1RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC98O2(i)=fC98O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C922O2 + PINALPA1RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC922O2(i)=fC922O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C923O2 + PINALPA1RO2 = C19H30O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C923O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC923O2(i)=fC923O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H30O7(i)=fC19H30O7(i)+1; 

i=i+1;
Rnames{i} = 'C924O2 + PINALPA1RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C924O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC924O2(i)=fC924O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C816CO3 + PINALPA1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C816CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC816CO3(i)=fC816CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'NORLIMO2 + PINALPA1RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NORLIMO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNORLIMO2(i)=fNORLIMO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'LMKAO2 + PINALPA1RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKAO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fLMKAO2(i)=fLMKAO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'LMKBO2 + PINALPA1RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKBO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fLMKBO2(i)=fLMKBO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C926O2 + PINALPA1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C926O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC926O2(i)=fC926O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C817CO3 + PINALPA1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC817CO3(i)=fC817CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'LMLKAO2 + PINALPA1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKAO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fLMLKAO2(i)=fLMLKAO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'LMLKBO2 + PINALPA1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMLKBO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fLMLKBO2(i)=fLMLKBO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C823CO3 + PINALPA1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC823CO3(i)=fC823CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C925O2 + PINALPA1RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C925O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC925O2(i)=fC925O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'NOPINAO2 + PINALPA1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINAO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNOPINAO2(i)=fNOPINAO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NOPINBO2 + PINALPA1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINBO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNOPINBO2(i)=fNOPINBO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NOPINCO2 + PINALPA1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINCO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNOPINCO2(i)=fNOPINCO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'NOPINDO2 + PINALPA1RO2 = C19H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINDO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNOPINDO2(i)=fNOPINDO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O7(i)=fC19H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C918O2 + PINALPA1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC918O2(i)=fC918O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C9DCO2 + PINALPA1RO2 = C19H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C9DCO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC9DCO2(i)=fC9DCO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H26O8(i)=fC19H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C915O2 + PINALPA1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C915O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC915O2(i)=fC915O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C917O2 + PINALPA1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C917O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC917O2(i)=fC917O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C919O2 + PINALPA1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C919O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC919O2(i)=fC919O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C914O2 + PINALPA1RO2 = C19H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C914O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC914O2(i)=fC914O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H26O9(i)=fC19H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C916O2 + PINALPA1RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C916O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC916O2(i)=fC916O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C88CO3 + PINALPA1RO2 = C19H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC88CO3(i)=fC88CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H26O9(i)=fC19H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C87CO3 + PINALPA1RO2 = C19H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C87CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC87CO3(i)=fC87CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H26O10(i)=fC19H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C822CO3 + PINALPA1RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C822CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC822CO3(i)=fC822CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'NLMKAO2 + PINALPA1RO2 = C19H29O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLMKAO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNLMKAO2(i)=fNLMKAO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC19H29O7NO3(i)=fC19H29O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'C85O2 + PINALPA1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC85O2(i)=fC85O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C89O2 + PINALPA1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC89O2(i)=fC89O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C86O2 + PINALPA1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC86O2(i)=fC86O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C811O2 + PINALPA1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC811O2(i)=fC811O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C810O2 + PINALPA1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC810O2(i)=fC810O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C812O2 + PINALPA1RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC812O2(i)=fC812O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C813O2 + PINALPA1RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC813O2(i)=fC813O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C729CO3 + PINALPA1RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C729CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC729CO3(i)=fC729CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C816O2 + PINALPA1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C816O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC816O2(i)=fC816O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C817O2 + PINALPA1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC817O2(i)=fC817O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C826O2 + PINALPA1RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C826O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC826O2(i)=fC826O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C822O2 + PINALPA1RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C822O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC822O2(i)=fC822O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C818O2 + PINALPA1RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C818O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC818O2(i)=fC818O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C823O2 + PINALPA1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC823O2(i)=fC823O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C819O2 + PINALPA1RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C819O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC819O2(i)=fC819O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C727CO3 + PINALPA1RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C727CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC727CO3(i)=fC727CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C731CO3 + PINALPA1RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C731CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC731CO3(i)=fC731CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C824O2 + PINALPA1RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C824O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC824O2(i)=fC824O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C820O2 + PINALPA1RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C820O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC820O2(i)=fC820O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C825O2 + PINALPA1RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C825O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC825O2(i)=fC825O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C732CO3 + PINALPA1RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C732CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC732CO3(i)=fC732CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C8BCO2 + PINALPA1RO2 = C18H28O6';
k(:,i) = kROORa.*fROOR;
Gstr{i,1} = 'C8BCO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC8BCO2(i)=fC8BCO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H28O6(i)=fC18H28O6(i)+1; 

i=i+1;
Rnames{i} = 'C88O2 + PINALPA1RO2 = C18H26O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC88O2(i)=fC88O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H26O8(i)=fC18H26O8(i)+1; 

i=i+1;
Rnames{i} = 'C718CO3 + PINALPA1RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C718CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC718CO3(i)=fC718CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C87O2 + PINALPA1RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C87O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fC87O2(i)=fC87O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C721CO3 + PINALPA1RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'PINALPA1RO2'; 
fC721CO3(i)=fC721CO3(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'NC826O2 + PINALPA1RO2 = C18H27O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC826O2'; Gstr{i,2} = 'PINALPA1RO2'; 
fNC826O2(i)=fNC826O2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fC18H27O8NO3(i)=fC18H27O8NO3(i)+1;

i=i+1;
Rnames{i} = 'BUT2OLO2 + PINALPA1RO2 = BUT2OLdimer';
k(:,i) = kROORb.*fROOR1;
Gstr{i,1} = 'BUT2OLO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fBUT2OLO2(i)=fBUT2OLO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1; 

i=i+1;
Rnames{i} = 'CHEXO2 + PINALPA1RO2 = CHEXdimer';
k(:,i) = kROORa.*fROOR2;
Gstr{i,1} = 'CHEXO2'; Gstr{i,2} = 'PINALPA1RO2'; 
fCHEXO2(i)=fCHEXO2(i)-1; fPINALPA1RO2(i)=fPINALPA1RO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1; 

%% APC10H15O7RO2 rxns
i=i+1;
Rnames{i} = 'APINAO2 + APC10H15O7RO2 = C20H32O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H32O8(i)=fC20H32O8(i)+1; 

i=i+1;
Rnames{i} = 'APINBO2 + APC10H15O7RO2 = C20H32O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H32O8(i)=fC20H32O8(i)+1; 

i=i+1;
Rnames{i} = 'APINCO2 + APC10H15O7RO2 = C20H32O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H32O8(i)=fC20H32O8(i)+1; 

i=i+1;
Rnames{i} = 'C107O2 + APC10H15O7RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC107O2(i)=fC107O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C109O2 + APC10H15O7RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC109O2(i)=fC109O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C106O2 + APC10H15O7RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC106O2(i)=fC106O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C920CO3 + APC10H15O7RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC920CO3(i)=fC920CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C108O2 + APC10H15O7RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC108O2(i)=fC108O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1; 

i=i+1;
Rnames{i} = 'PINALO2 + APC10H15O7RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fPINALO2(i)=fPINALO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C96CO3 + APC10H15O7RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC96CO3(i)=fC96CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C923CO3 + APC10H15O7RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C923CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC923CO3(i)=fC923CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'LIMAO2 + APC10H15O7RO2 = C20H32O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMAO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fLIMAO2(i)=fLIMAO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H32O8(i)=fC20H32O8(i)+1; 

i=i+1;
Rnames{i} = 'LIMALBO2 + APC10H15O7RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALBO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fLIMALBO2(i)=fLIMALBO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'LIMCO2 + APC10H15O7RO2 = C20H32O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMCO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fLIMCO2(i)=fLIMCO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H32O8(i)=fC20H32O8(i)+1; 

i=i+1;
Rnames{i} = 'LIMALO2 + APC10H15O7RO2 = C20H32O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMALO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fLIMALO2(i)=fLIMALO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H32O10(i)=fC20H32O10(i)+1; 

i=i+1;
Rnames{i} = 'LIMBO2 + APC10H15O7RO2 = C20H32O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMBO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fLIMBO2(i)=fLIMBO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H32O8(i)=fC20H32O8(i)+1; 

i=i+1;
Rnames{i} = 'LIMALAO2 + APC10H15O7RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMALAO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fLIMALAO2(i)=fLIMALAO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'NAPINAO2 + APC10H15O7RO2 = C20H31O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NAPINAO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNAPINAO2(i)=fNAPINAO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H31O7NO3(i)=fC20H31O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NAPINBO2 + APC10H15O7RO2 = C20H31O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NAPINBO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNAPINBO2(i)=fNAPINBO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H31O7NO3(i)=fC20H31O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'BPINAO2 + APC10H15O7RO2 = C20H32O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINAO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fBPINAO2(i)=fBPINAO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H32O8(i)=fC20H32O8(i)+1; 

i=i+1;
Rnames{i} = 'BPINBO2 + APC10H15O7RO2 = C20H32O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINBO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fBPINBO2(i)=fBPINBO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H32O8(i)=fC20H32O8(i)+1; 

i=i+1;
Rnames{i} = 'BPINCO2 + APC10H15O7RO2 = C20H32O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINCO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fBPINCO2(i)=fBPINCO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H32O8(i)=fC20H32O8(i)+1; 

i=i+1;
Rnames{i} = 'C918CO3 + APC10H15O7RO2 = C20H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC918CO3(i)=fC918CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H30O9(i)=fC20H30O9(i)+1; 

i=i+1;
Rnames{i} = 'NC102O2 + APC10H15O7RO2 = C20H29O9NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC102O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNC102O2(i)=fNC102O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H29O9NO3(i)=fC20H29O9NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC101O2 + APC10H15O7RO2 = C20H29O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC101O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNC101O2(i)=fNC101O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H29O8NO3(i)=fC20H29O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMO2 + APC10H15O7RO2 = C20H31O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLIMO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNLIMO2(i)=fNLIMO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H31O7NO3(i)=fC20H31O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMALO2 + APC10H15O7RO2 = C20H31O9NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLIMALO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNLIMALO2(i)=fNLIMALO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H31O9NO3(i)=fC20H31O9NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC91CO3 + APC10H15O7RO2 = C20H29O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC91CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNC91CO3(i)=fNC91CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H29O8NO3(i)=fC20H29O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINAO2 + APC10H15O7RO2 = C20H31O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NBPINAO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNBPINAO2(i)=fNBPINAO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H31O7NO3(i)=fC20H31O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINBO2 + APC10H15O7RO2 = C20H31O7NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NBPINBO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNBPINBO2(i)=fNBPINBO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC20H31O7NO3(i)=fC20H31O7NO3(i)+1; 

i=i+1;
Rnames{i} = 'C96O2 + APC10H15O7RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC96O2(i)=fC96O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C89CO3 + APC10H15O7RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC89CO3(i)=fC89CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C920O2 + APC10H15O7RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC920O2(i)=fC920O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C97O2 + APC10H15O7RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC97O2(i)=fC97O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C85CO3 + APC10H15O7RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC85CO3(i)=fC85CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C811CO3 + APC10H15O7RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC811CO3(i)=fC811CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C921O2 + APC10H15O7RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC921O2(i)=fC921O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C98O2 + APC10H15O7RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC98O2(i)=fC98O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C922O2 + APC10H15O7RO2 = C19H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC922O2(i)=fC922O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H30O11(i)=fC19H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C923O2 + APC10H15O7RO2 = C19H30O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C923O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC923O2(i)=fC923O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H30O8(i)=fC19H30O8(i)+1; 

i=i+1;
Rnames{i} = 'C924O2 + APC10H15O7RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C924O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC924O2(i)=fC924O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C816CO3 + APC10H15O7RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C816CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC816CO3(i)=fC816CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'NORLIMO2 + APC10H15O7RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NORLIMO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNORLIMO2(i)=fNORLIMO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'LMKAO2 + APC10H15O7RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKAO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fLMKAO2(i)=fLMKAO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'LMKBO2 + APC10H15O7RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LMKBO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fLMKBO2(i)=fLMKBO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C926O2 + APC10H15O7RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C926O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC926O2(i)=fC926O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C817CO3 + APC10H15O7RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C817CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC817CO3(i)=fC817CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'LMLKAO2 + APC10H15O7RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LMLKAO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fLMLKAO2(i)=fLMLKAO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'LMLKBO2 + APC10H15O7RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LMLKBO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fLMLKBO2(i)=fLMLKBO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C823CO3 + APC10H15O7RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C823CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC823CO3(i)=fC823CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C925O2 + APC10H15O7RO2 = C19H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C925O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC925O2(i)=fC925O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H30O11(i)=fC19H30O11(i)+1; 

i=i+1;
Rnames{i} = 'NOPINAO2 + APC10H15O7RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINAO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNOPINAO2(i)=fNOPINAO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'NOPINBO2 + APC10H15O7RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINBO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNOPINBO2(i)=fNOPINBO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'NOPINCO2 + APC10H15O7RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINCO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNOPINCO2(i)=fNOPINCO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'NOPINDO2 + APC10H15O7RO2 = C19H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINDO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNOPINDO2(i)=fNOPINDO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O8(i)=fC19H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C918O2 + APC10H15O7RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C918O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC918O2(i)=fC918O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C9DCO2 + APC10H15O7RO2 = C19H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C9DCO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC9DCO2(i)=fC9DCO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H26O9(i)=fC19H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C915O2 + APC10H15O7RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C915O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC915O2(i)=fC915O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C917O2 + APC10H15O7RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C917O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC917O2(i)=fC917O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C919O2 + APC10H15O7RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C919O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC919O2(i)=fC919O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C914O2 + APC10H15O7RO2 = C19H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C914O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC914O2(i)=fC914O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H26O10(i)=fC19H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C916O2 + APC10H15O7RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C916O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC916O2(i)=fC916O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C88CO3 + APC10H15O7RO2 = C19H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C88CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC88CO3(i)=fC88CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H26O10(i)=fC19H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C87CO3 + APC10H15O7RO2 = C19H26O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C87CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC87CO3(i)=fC87CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H26O11(i)=fC19H26O11(i)+1; 

i=i+1;
Rnames{i} = 'C822CO3 + APC10H15O7RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C822CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC822CO3(i)=fC822CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'NLMKAO2 + APC10H15O7RO2 = C19H29O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLMKAO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNLMKAO2(i)=fNLMKAO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC19H29O8NO3(i)=fC19H29O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'C85O2 + APC10H15O7RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC85O2(i)=fC85O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C89O2 + APC10H15O7RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC89O2(i)=fC89O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C86O2 + APC10H15O7RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC86O2(i)=fC86O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C811O2 + APC10H15O7RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC811O2(i)=fC811O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C810O2 + APC10H15O7RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC810O2(i)=fC810O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C812O2 + APC10H15O7RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC812O2(i)=fC812O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C813O2 + APC10H15O7RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC813O2(i)=fC813O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C729CO3 + APC10H15O7RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C729CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC729CO3(i)=fC729CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C816O2 + APC10H15O7RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C816O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC816O2(i)=fC816O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C817O2 + APC10H15O7RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C817O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC817O2(i)=fC817O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C826O2 + APC10H15O7RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C826O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC826O2(i)=fC826O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C822O2 + APC10H15O7RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C822O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC822O2(i)=fC822O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C818O2 + APC10H15O7RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C818O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC818O2(i)=fC818O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C823O2 + APC10H15O7RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C823O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC823O2(i)=fC823O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C819O2 + APC10H15O7RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C819O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC819O2(i)=fC819O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C727CO3 + APC10H15O7RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C727CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC727CO3(i)=fC727CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C731CO3 + APC10H15O7RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C731CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC731CO3(i)=fC731CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C824O2 + APC10H15O7RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C824O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC824O2(i)=fC824O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C820O2 + APC10H15O7RO2 = C18H26O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C820O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC820O2(i)=fC820O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H26O11(i)=fC18H26O11(i)+1; 

i=i+1;
Rnames{i} = 'C825O2 + APC10H15O7RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C825O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC825O2(i)=fC825O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C732CO3 + APC10H15O7RO2 = C18H26O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C732CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC732CO3(i)=fC732CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H26O11(i)=fC18H26O11(i)+1; 

i=i+1;
Rnames{i} = 'C8BCO2 + APC10H15O7RO2 = C18H28O7';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C8BCO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC8BCO2(i)=fC8BCO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H28O7(i)=fC18H28O7(i)+1; 

i=i+1;
Rnames{i} = 'C88O2 + APC10H15O7RO2 = C18H26O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C88O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC88O2(i)=fC88O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H26O9(i)=fC18H26O9(i)+1; 

i=i+1;
Rnames{i} = 'C718CO3 + APC10H15O7RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C718CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC718CO3(i)=fC718CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C87O2 + APC10H15O7RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C87O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC87O2(i)=fC87O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C721CO3 + APC10H15O7RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'APC10H15O7RO2'; 
fC721CO3(i)=fC721CO3(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'NC826O2 + APC10H15O7RO2 = C18H27O9NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC826O2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fNC826O2(i)=fNC826O2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fC18H27O9NO3(i)=fC18H27O9NO3(i)+1; 

i=i+1;
Rnames{i} = 'BUT2OLO2 + APC10H15O7RO2 = BUT2OLdimer';
k(:,i) = kROORb.*fROOR1;
Gstr{i,1} = 'BUT2OLO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fBUT2OLO2(i)=fBUT2OLO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1; 

i=i+1;
Rnames{i} = 'CHEXO2 + APC10H15O7RO2 = CHEXdimer';
k(:,i) = kROORb.*fROOR2;
Gstr{i,1} = 'CHEXO2'; Gstr{i,2} = 'APC10H15O7RO2'; 
fCHEXO2(i)=fCHEXO2(i)-1; fAPC10H15O7RO2(i)=fAPC10H15O7RO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1; 

%% APC10H15O8RO2 rxns
i=i+1;
Rnames{i} = 'APINAO2 + APC10H15O8RO2 = C20H32O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H32O9(i)=fC20H32O9(i)+1; 

i=i+1;
Rnames{i} = 'APINBO2 + APC10H15O8RO2 = C20H32O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H32O9(i)=fC20H32O9(i)+1; 

i=i+1;
Rnames{i} = 'APINCO2 + APC10H15O8RO2 = C20H32O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H32O9(i)=fC20H32O9(i)+1; 

i=i+1;
Rnames{i} = 'C107O2 + APC10H15O8RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC107O2(i)=fC107O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C109O2 + APC10H15O8RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC109O2(i)=fC109O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C106O2 + APC10H15O8RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC106O2(i)=fC106O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C920CO3 + APC10H15O8RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC920CO3(i)=fC920CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C108O2 + APC10H15O8RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC108O2(i)=fC108O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1; 

i=i+1;
Rnames{i} = 'PINALO2 + APC10H15O8RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fPINALO2(i)=fPINALO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C96CO3 + APC10H15O8RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC96CO3(i)=fC96CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C923CO3 + APC10H15O8RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C923CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC923CO3(i)=fC923CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1; 

i=i+1;
Rnames{i} = 'LIMAO2 + APC10H15O8RO2 = C20H32O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMAO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fLIMAO2(i)=fLIMAO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H32O9(i)=fC20H32O9(i)+1; 

i=i+1;
Rnames{i} = 'LIMALBO2 + APC10H15O8RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMALBO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fLIMALBO2(i)=fLIMALBO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1; 

i=i+1;
Rnames{i} = 'LIMCO2 + APC10H15O8RO2 = C20H32O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMCO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fLIMCO2(i)=fLIMCO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H32O9(i)=fC20H32O9(i)+1; 

i=i+1;
Rnames{i} = 'LIMALO2 + APC10H15O8RO2 = C20H32O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMALO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fLIMALO2(i)=fLIMALO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H32O11(i)=fC20H32O11(i)+1; 

i=i+1;
Rnames{i} = 'LIMBO2 + APC10H15O8RO2 = C20H32O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'LIMBO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fLIMBO2(i)=fLIMBO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H32O9(i)=fC20H32O9(i)+1; 

i=i+1;
Rnames{i} = 'LIMALAO2 + APC10H15O8RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMALAO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fLIMALAO2(i)=fLIMALAO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1; 

i=i+1;
Rnames{i} = 'NAPINAO2 + APC10H15O8RO2 = C20H31O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NAPINAO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNAPINAO2(i)=fNAPINAO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H31O8NO3(i)=fC20H31O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'NAPINBO2 + APC10H15O8RO2 = C20H31O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NAPINBO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNAPINBO2(i)=fNAPINBO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H31O8NO3(i)=fC20H31O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'BPINAO2 + APC10H15O8RO2 = C20H32O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINAO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fBPINAO2(i)=fBPINAO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H32O9(i)=fC20H32O9(i)+1; 

i=i+1;
Rnames{i} = 'BPINBO2 + APC10H15O8RO2 = C20H32O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINBO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fBPINBO2(i)=fBPINBO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H32O9(i)=fC20H32O9(i)+1; 

i=i+1;
Rnames{i} = 'BPINCO2 + APC10H15O8RO2 = C20H32O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'BPINCO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fBPINCO2(i)=fBPINCO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H32O9(i)=fC20H32O9(i)+1; 

i=i+1;
Rnames{i} = 'C918CO3 + APC10H15O8RO2 = C20H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C918CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC918CO3(i)=fC918CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H30O10(i)=fC20H30O10(i)+1; 

i=i+1;
Rnames{i} = 'NC102O2 + APC10H15O8RO2 = C20H29O10NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NC102O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNC102O2(i)=fNC102O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H29O10NO3(i)=fC20H29O10NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC101O2 + APC10H15O8RO2 = C20H29O9NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC101O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNC101O2(i)=fNC101O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H29O9NO3(i)=fC20H29O9NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMO2 + APC10H15O8RO2 = C20H31O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLIMO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNLIMO2(i)=fNLIMO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H31O8NO3(i)=fC20H31O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMALO2 + APC10H15O8RO2 = C20H31O10NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NLIMALO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNLIMALO2(i)=fNLIMALO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H31O10NO3(i)=fC20H31O10NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC91CO3 + APC10H15O8RO2 = C20H29O9NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NC91CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNC91CO3(i)=fNC91CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H29O9NO3(i)=fC20H29O9NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINAO2 + APC10H15O8RO2 = C20H31O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NBPINAO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNBPINAO2(i)=fNBPINAO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H31O8NO3(i)=fC20H31O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINBO2 + APC10H15O8RO2 = C20H31O8NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NBPINBO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNBPINBO2(i)=fNBPINBO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC20H31O8NO3(i)=fC20H31O8NO3(i)+1; 

i=i+1;
Rnames{i} = 'C96O2 + APC10H15O8RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC96O2(i)=fC96O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C89CO3 + APC10H15O8RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC89CO3(i)=fC89CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C920O2 + APC10H15O8RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC920O2(i)=fC920O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C97O2 + APC10H15O8RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC97O2(i)=fC97O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C85CO3 + APC10H15O8RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC85CO3(i)=fC85CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C811CO3 + APC10H15O8RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC811CO3(i)=fC811CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C921O2 + APC10H15O8RO2 = C19H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC921O2(i)=fC921O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H30O11(i)=fC19H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C98O2 + APC10H15O8RO2 = C19H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC98O2(i)=fC98O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H30O11(i)=fC19H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C922O2 + APC10H15O8RO2 = C19H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC922O2(i)=fC922O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H30O12(i)=fC19H30O12(i)+1; 

i=i+1;
Rnames{i} = 'C923O2 + APC10H15O8RO2 = C19H30O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C923O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC923O2(i)=fC923O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H30O9(i)=fC19H30O9(i)+1; 

i=i+1;
Rnames{i} = 'C924O2 + APC10H15O8RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C924O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC924O2(i)=fC924O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C816CO3 + APC10H15O8RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C816CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC816CO3(i)=fC816CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'NORLIMO2 + APC10H15O8RO2 = C19H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NORLIMO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNORLIMO2(i)=fNORLIMO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H30O11(i)=fC19H30O11(i)+1; 

i=i+1;
Rnames{i} = 'LMKAO2 + APC10H15O8RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LMKAO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fLMKAO2(i)=fLMKAO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'LMKBO2 + APC10H15O8RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LMKBO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fLMKBO2(i)=fLMKBO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C926O2 + APC10H15O8RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C926O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC926O2(i)=fC926O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C817CO3 + APC10H15O8RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C817CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC817CO3(i)=fC817CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'LMLKAO2 + APC10H15O8RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LMLKAO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fLMLKAO2(i)=fLMLKAO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'LMLKBO2 + APC10H15O8RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LMLKBO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fLMLKBO2(i)=fLMLKBO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C823CO3 + APC10H15O8RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C823CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC823CO3(i)=fC823CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C925O2 + APC10H15O8RO2 = C19H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C925O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC925O2(i)=fC925O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H30O12(i)=fC19H30O12(i)+1; 

i=i+1;
Rnames{i} = 'NOPINAO2 + APC10H15O8RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINAO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNOPINAO2(i)=fNOPINAO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'NOPINBO2 + APC10H15O8RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINBO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNOPINBO2(i)=fNOPINBO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'NOPINCO2 + APC10H15O8RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINCO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNOPINCO2(i)=fNOPINCO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'NOPINDO2 + APC10H15O8RO2 = C19H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'NOPINDO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNOPINDO2(i)=fNOPINDO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O9(i)=fC19H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C918O2 + APC10H15O8RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C918O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC918O2(i)=fC918O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C9DCO2 + APC10H15O8RO2 = C19H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C9DCO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC9DCO2(i)=fC9DCO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H26O10(i)=fC19H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C915O2 + APC10H15O8RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C915O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC915O2(i)=fC915O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C917O2 + APC10H15O8RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C917O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC917O2(i)=fC917O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C919O2 + APC10H15O8RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C919O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC919O2(i)=fC919O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C914O2 + APC10H15O8RO2 = C19H26O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C914O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC914O2(i)=fC914O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H26O11(i)=fC19H26O11(i)+1; 

i=i+1;
Rnames{i} = 'C916O2 + APC10H15O8RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C916O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC916O2(i)=fC916O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C88CO3 + APC10H15O8RO2 = C19H26O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C88CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC88CO3(i)=fC88CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H26O11(i)=fC19H26O11(i)+1; 

i=i+1;
Rnames{i} = 'C87CO3 + APC10H15O8RO2 = C19H26O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C87CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC87CO3(i)=fC87CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H26O12(i)=fC19H26O12(i)+1; 

i=i+1;
Rnames{i} = 'C822CO3 + APC10H15O8RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C822CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC822CO3(i)=fC822CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'NLMKAO2 + APC10H15O8RO2 = C19H29O9NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLMKAO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNLMKAO2(i)=fNLMKAO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC19H29O9NO3(i)=fC19H29O9NO3(i)+1; 

i=i+1;
Rnames{i} = 'C85O2 + APC10H15O8RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC85O2(i)=fC85O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C89O2 + APC10H15O8RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC89O2(i)=fC89O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C86O2 + APC10H15O8RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC86O2(i)=fC86O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C811O2 + APC10H15O8RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC811O2(i)=fC811O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C810O2 + APC10H15O8RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC810O2(i)=fC810O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C812O2 + APC10H15O8RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC812O2(i)=fC812O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C813O2 + APC10H15O8RO2 = C18H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC813O2(i)=fC813O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O12(i)=fC18H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C729CO3 + APC10H15O8RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C729CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC729CO3(i)=fC729CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C816O2 + APC10H15O8RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C816O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC816O2(i)=fC816O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C817O2 + APC10H15O8RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C817O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC817O2(i)=fC817O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C826O2 + APC10H15O8RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C826O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC826O2(i)=fC826O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C822O2 + APC10H15O8RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C822O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC822O2(i)=fC822O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C818O2 + APC10H15O8RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C818O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC818O2(i)=fC818O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C823O2 + APC10H15O8RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C823O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC823O2(i)=fC823O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C819O2 + APC10H15O8RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C819O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC819O2(i)=fC819O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C727CO3 + APC10H15O8RO2 = C18H26O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C727CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC727CO3(i)=fC727CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H26O11(i)=fC18H26O11(i)+1; 

i=i+1;
Rnames{i} = 'C731CO3 + APC10H15O8RO2 = C18H26O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C731CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC731CO3(i)=fC731CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H26O11(i)=fC18H26O11(i)+1; 

i=i+1;
Rnames{i} = 'C824O2 + APC10H15O8RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C824O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC824O2(i)=fC824O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C820O2 + APC10H15O8RO2 = C18H26O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C820O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC820O2(i)=fC820O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H26O12(i)=fC18H26O12(i)+1; 

i=i+1;
Rnames{i} = 'C825O2 + APC10H15O8RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C825O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC825O2(i)=fC825O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C732CO3 + APC10H15O8RO2 = C18H26O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C732CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC732CO3(i)=fC732CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H26O12(i)=fC18H26O12(i)+1; 

i=i+1;
Rnames{i} = 'C8BCO2 + APC10H15O8RO2 = C18H28O8';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C8BCO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC8BCO2(i)=fC8BCO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H28O8(i)=fC18H28O8(i)+1; 

i=i+1;
Rnames{i} = 'C88O2 + APC10H15O8RO2 = C18H26O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C88O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC88O2(i)=fC88O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H26O10(i)=fC18H26O10(i)+1; 

i=i+1;
Rnames{i} = 'C718CO3 + APC10H15O8RO2 = C18H26O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C718CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC718CO3(i)=fC718CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H26O11(i)=fC18H26O11(i)+1; 

i=i+1;
Rnames{i} = 'C87O2 + APC10H15O8RO2 = C18H26O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C87O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC87O2(i)=fC87O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H26O11(i)=fC18H26O11(i)+1; 

i=i+1;
Rnames{i} = 'C721CO3 + APC10H15O8RO2 = C18H26O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'APC10H15O8RO2'; 
fC721CO3(i)=fC721CO3(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H26O11(i)=fC18H26O11(i)+1; 

i=i+1;
Rnames{i} = 'NC826O2 + APC10H15O8RO2 = C18H27O10NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NC826O2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fNC826O2(i)=fNC826O2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fC18H27O10NO3(i)=fC18H27O10NO3(i)+1; 

i=i+1;
Rnames{i} = 'BUT2OLO2 + APC10H15O8RO2 = BUT2OLdimer';
k(:,i) = kROORb.*fROOR1;
Gstr{i,1} = 'BUT2OLO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fBUT2OLO2(i)=fBUT2OLO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1; 

i=i+1;
Rnames{i} = 'CHEXO2 + APC10H15O8RO2 = CHEXdimer';
k(:,i) = kROORb.*fROOR2;
Gstr{i,1} = 'CHEXO2'; Gstr{i,2} = 'APC10H15O8RO2'; 
fCHEXO2(i)=fCHEXO2(i)-1; fAPC10H15O8RO2(i)=fAPC10H15O8RO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1; 

%% APC10H15O9RO2 rxns
i=i+1;
Rnames{i} = 'APINAO2 + APC10H15O9RO2 = C20H32O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H32O10(i)=fC20H32O10(i)+1; 

i=i+1;
Rnames{i} = 'APINBO2 + APC10H15O9RO2 = C20H32O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H32O10(i)=fC20H32O10(i)+1; 

i=i+1;
Rnames{i} = 'APINCO2 + APC10H15O9RO2 = C20H32O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H32O10(i)=fC20H32O10(i)+1; 

i=i+1;
Rnames{i} = 'C107O2 + APC10H15O9RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC107O2(i)=fC107O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C109O2 + APC10H15O9RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC109O2(i)=fC109O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C106O2 + APC10H15O9RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC106O2(i)=fC106O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1; 

i=i+1;
Rnames{i} = 'C920CO3 + APC10H15O9RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC920CO3(i)=fC920CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1; 

i=i+1;
Rnames{i} = 'C108O2 + APC10H15O9RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC108O2(i)=fC108O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1; 

i=i+1;
Rnames{i} = 'PINALO2 + APC10H15O9RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fPINALO2(i)=fPINALO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C96CO3 + APC10H15O9RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC96CO3(i)=fC96CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C923CO3 + APC10H15O9RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C923CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC923CO3(i)=fC923CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1; 

i=i+1;
Rnames{i} = 'LIMAO2 + APC10H15O9RO2 = C20H32O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMAO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fLIMAO2(i)=fLIMAO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H32O10(i)=fC20H32O10(i)+1; 

i=i+1;
Rnames{i} = 'LIMALBO2 + APC10H15O9RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMALBO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fLIMALBO2(i)=fLIMALBO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1; 

i=i+1;
Rnames{i} = 'LIMCO2 + APC10H15O9RO2 = C20H32O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMCO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fLIMCO2(i)=fLIMCO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H32O10(i)=fC20H32O10(i)+1; 

i=i+1;
Rnames{i} = 'LIMALO2 + APC10H15O9RO2 = C20H32O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMALO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fLIMALO2(i)=fLIMALO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H32O12(i)=fC20H32O12(i)+1; 

i=i+1;
Rnames{i} = 'LIMBO2 + APC10H15O9RO2 = C20H32O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMBO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fLIMBO2(i)=fLIMBO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H32O10(i)=fC20H32O10(i)+1; 

i=i+1;
Rnames{i} = 'LIMALAO2 + APC10H15O9RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMALAO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fLIMALAO2(i)=fLIMALAO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1; 

i=i+1;
Rnames{i} = 'NAPINAO2 + APC10H15O9RO2 = C20H31O9NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NAPINAO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNAPINAO2(i)=fNAPINAO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H31O9NO3(i)=fC20H31O9NO3(i)+1; 

i=i+1;
Rnames{i} = 'NAPINBO2 + APC10H15O9RO2 = C20H31O9NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NAPINBO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNAPINBO2(i)=fNAPINBO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H31O9NO3(i)=fC20H31O9NO3(i)+1; 

i=i+1;
Rnames{i} = 'BPINAO2 + APC10H15O9RO2 = C20H32O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'BPINAO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fBPINAO2(i)=fBPINAO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H32O10(i)=fC20H32O10(i)+1; 

i=i+1;
Rnames{i} = 'BPINBO2 + APC10H15O9RO2 = C20H32O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'BPINBO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fBPINBO2(i)=fBPINBO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H32O10(i)=fC20H32O10(i)+1; 

i=i+1;
Rnames{i} = 'BPINCO2 + APC10H15O9RO2 = C20H32O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'BPINCO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fBPINCO2(i)=fBPINCO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H32O10(i)=fC20H32O10(i)+1; 

i=i+1;
Rnames{i} = 'C918CO3 + APC10H15O9RO2 = C20H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C918CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC918CO3(i)=fC918CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H30O11(i)=fC20H30O11(i)+1; 

i=i+1;
Rnames{i} = 'NC102O2 + APC10H15O9RO2 = C20H29O11NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NC102O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNC102O2(i)=fNC102O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H29O11NO3(i)=fC20H29O11NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC101O2 + APC10H15O9RO2 = C20H29O10NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NC101O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNC101O2(i)=fNC101O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H29O10NO3(i)=fC20H29O10NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMO2 + APC10H15O9RO2 = C20H31O9NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NLIMO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNLIMO2(i)=fNLIMO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H31O9NO3(i)=fC20H31O9NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMALO2 + APC10H15O9RO2 = C20H31O11NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NLIMALO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNLIMALO2(i)=fNLIMALO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H31O11NO3(i)=fC20H31O11NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC91CO3 + APC10H15O9RO2 = C20H29O10NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NC91CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNC91CO3(i)=fNC91CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H29O10NO3(i)=fC20H29O10NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINAO2 + APC10H15O9RO2 = C20H31O9NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NBPINAO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNBPINAO2(i)=fNBPINAO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H31O9NO3(i)=fC20H31O9NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINBO2 + APC10H15O9RO2 = C20H31O9NO3';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NBPINBO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNBPINBO2(i)=fNBPINBO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC20H31O9NO3(i)=fC20H31O9NO3(i)+1; 

i=i+1;
Rnames{i} = 'C96O2 + APC10H15O9RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC96O2(i)=fC96O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C89CO3 + APC10H15O9RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC89CO3(i)=fC89CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C920O2 + APC10H15O9RO2 = C19H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC920O2(i)=fC920O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H30O11(i)=fC19H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C97O2 + APC10H15O9RO2 = C19H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC97O2(i)=fC97O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H30O11(i)=fC19H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C85CO3 + APC10H15O9RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC85CO3(i)=fC85CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C811CO3 + APC10H15O9RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC811CO3(i)=fC811CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C921O2 + APC10H15O9RO2 = C19H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC921O2(i)=fC921O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H30O12(i)=fC19H30O12(i)+1; 

i=i+1;
Rnames{i} = 'C98O2 + APC10H15O9RO2 = C19H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC98O2(i)=fC98O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H30O12(i)=fC19H30O12(i)+1; 

i=i+1;
Rnames{i} = 'C922O2 + APC10H15O9RO2 = C19H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC922O2(i)=fC922O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H30O13(i)=fC19H30O13(i)+1; 

i=i+1;
Rnames{i} = 'C923O2 + APC10H15O9RO2 = C19H30O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C923O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC923O2(i)=fC923O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H30O10(i)=fC19H30O10(i)+1; 

i=i+1;
Rnames{i} = 'C924O2 + APC10H15O9RO2 = C19H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C924O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC924O2(i)=fC924O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H30O11(i)=fC19H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C816CO3 + APC10H15O9RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C816CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC816CO3(i)=fC816CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'NORLIMO2 + APC10H15O9RO2 = C19H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NORLIMO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNORLIMO2(i)=fNORLIMO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H30O12(i)=fC19H30O12(i)+1; 

i=i+1;
Rnames{i} = 'LMKAO2 + APC10H15O9RO2 = C19H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LMKAO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fLMKAO2(i)=fLMKAO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H30O11(i)=fC19H30O11(i)+1; 

i=i+1;
Rnames{i} = 'LMKBO2 + APC10H15O9RO2 = C19H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LMKBO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fLMKBO2(i)=fLMKBO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H30O11(i)=fC19H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C926O2 + APC10H15O9RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C926O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC926O2(i)=fC926O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C817CO3 + APC10H15O9RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C817CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC817CO3(i)=fC817CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'LMLKAO2 + APC10H15O9RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LMLKAO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fLMLKAO2(i)=fLMLKAO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'LMLKBO2 + APC10H15O9RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LMLKBO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fLMLKBO2(i)=fLMLKBO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C823CO3 + APC10H15O9RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C823CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC823CO3(i)=fC823CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C925O2 + APC10H15O9RO2 = C19H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C925O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC925O2(i)=fC925O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H30O13(i)=fC19H30O13(i)+1; 

i=i+1;
Rnames{i} = 'NOPINAO2 + APC10H15O9RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NOPINAO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNOPINAO2(i)=fNOPINAO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'NOPINBO2 + APC10H15O9RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NOPINBO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNOPINBO2(i)=fNOPINBO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'NOPINCO2 + APC10H15O9RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NOPINCO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNOPINCO2(i)=fNOPINCO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'NOPINDO2 + APC10H15O9RO2 = C19H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NOPINDO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNOPINDO2(i)=fNOPINDO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O10(i)=fC19H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C918O2 + APC10H15O9RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C918O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC918O2(i)=fC918O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C9DCO2 + APC10H15O9RO2 = C19H26O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C9DCO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC9DCO2(i)=fC9DCO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H26O11(i)=fC19H26O11(i)+1; 

i=i+1;
Rnames{i} = 'C915O2 + APC10H15O9RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C915O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC915O2(i)=fC915O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C917O2 + APC10H15O9RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C917O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC917O2(i)=fC917O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C919O2 + APC10H15O9RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C919O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC919O2(i)=fC919O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C914O2 + APC10H15O9RO2 = C19H26O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C914O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC914O2(i)=fC914O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H26O12(i)=fC19H26O12(i)+1; 

i=i+1;
Rnames{i} = 'C916O2 + APC10H15O9RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C916O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC916O2(i)=fC916O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C88CO3 + APC10H15O9RO2 = C19H26O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C88CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC88CO3(i)=fC88CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H26O12(i)=fC19H26O12(i)+1; 

i=i+1;
Rnames{i} = 'C87CO3 + APC10H15O9RO2 = C19H26O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C87CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC87CO3(i)=fC87CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H26O13(i)=fC19H26O13(i)+1;

i=i+1;
Rnames{i} = 'C822CO3 + APC10H15O9RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C822CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC822CO3(i)=fC822CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'NLMKAO2 + APC10H15O9RO2 = C19H29O10NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NLMKAO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNLMKAO2(i)=fNLMKAO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC19H29O10NO3(i)=fC19H29O10NO3(i)+1; 

i=i+1;
Rnames{i} = 'C85O2 + APC10H15O9RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC85O2(i)=fC85O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C89O2 + APC10H15O9RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC89O2(i)=fC89O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C86O2 + APC10H15O9RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC86O2(i)=fC86O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C811O2 + APC10H15O9RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC811O2(i)=fC811O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C810O2 + APC10H15O9RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC810O2(i)=fC810O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C812O2 + APC10H15O9RO2 = C18H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC812O2(i)=fC812O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O12(i)=fC18H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C813O2 + APC10H15O9RO2 = C18H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC813O2(i)=fC813O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O13(i)=fC18H28O13(i)+1; 

i=i+1;
Rnames{i} = 'C729CO3 + APC10H15O9RO2 = C18H26O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C729CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC729CO3(i)=fC729CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H26O11(i)=fC18H26O11(i)+1; 

i=i+1;
Rnames{i} = 'C816O2 + APC10H15O9RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C816O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC816O2(i)=fC816O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C817O2 + APC10H15O9RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C817O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC817O2(i)=fC817O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C826O2 + APC10H15O9RO2 = C18H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C826O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC826O2(i)=fC826O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O12(i)=fC18H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C822O2 + APC10H15O9RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C822O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC822O2(i)=fC822O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C818O2 + APC10H15O9RO2 = C18H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C818O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC818O2(i)=fC818O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O12(i)=fC18H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C823O2 + APC10H15O9RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C823O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC823O2(i)=fC823O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C819O2 + APC10H15O9RO2 = C18H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C819O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC819O2(i)=fC819O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O12(i)=fC18H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C727CO3 + APC10H15O9RO2 = C18H26O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C727CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC727CO3(i)=fC727CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H26O12(i)=fC18H26O12(i)+1; 

i=i+1;
Rnames{i} = 'C731CO3 + APC10H15O9RO2 = C18H26O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C731CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC731CO3(i)=fC731CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H26O12(i)=fC18H26O12(i)+1; 

i=i+1;
Rnames{i} = 'C824O2 + APC10H15O9RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C824O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC824O2(i)=fC824O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C820O2 + APC10H15O9RO2 = C18H26O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C820O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC820O2(i)=fC820O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H26O13(i)=fC18H26O13(i)+1; 

i=i+1;
Rnames{i} = 'C825O2 + APC10H15O9RO2 = C18H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C825O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC825O2(i)=fC825O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O12(i)=fC18H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C732CO3 + APC10H15O9RO2 = C18H26O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C732CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC732CO3(i)=fC732CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H26O13(i)=fC18H26O13(i)+1; 

i=i+1;
Rnames{i} = 'C8BCO2 + APC10H15O9RO2 = C18H28O9';
k(:,i) = kROORb.*fROOR;
Gstr{i,1} = 'C8BCO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC8BCO2(i)=fC8BCO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H28O9(i)=fC18H28O9(i)+1; 

i=i+1;
Rnames{i} = 'C88O2 + APC10H15O9RO2 = C18H26O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C88O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC88O2(i)=fC88O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H26O11(i)=fC18H26O11(i)+1; 

i=i+1;
Rnames{i} = 'C718CO3 + APC10H15O9RO2 = C18H26O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C718CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC718CO3(i)=fC718CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H26O12(i)=fC18H26O12(i)+1; 

i=i+1;
Rnames{i} = 'C87O2 + APC10H15O9RO2 = C18H26O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C87O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC87O2(i)=fC87O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H26O12(i)=fC18H26O12(i)+1; 

i=i+1;
Rnames{i} = 'C721CO3 + APC10H15O9RO2 = C18H26O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'APC10H15O9RO2'; 
fC721CO3(i)=fC721CO3(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H26O12(i)=fC18H26O12(i)+1; 

i=i+1;
Rnames{i} = 'NC826O2 + APC10H15O9RO2 = C18H27O11NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NC826O2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fNC826O2(i)=fNC826O2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fC18H27O11NO3(i)=fC18H27O11NO3(i)+1; 

i=i+1;
Rnames{i} = 'BUT2OLO2 + APC10H15O9RO2 = BUT2OLdimer';
k(:,i) = kROORc.*fROOR1;
Gstr{i,1} = 'BUT2OLO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fBUT2OLO2(i)=fBUT2OLO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1; 

i=i+1;
Rnames{i} = 'CHEXO2 + APC10H15O9RO2 = CHEXdimer';
k(:,i) = kROORb.*fROOR2;
Gstr{i,1} = 'CHEXO2'; Gstr{i,2} = 'APC10H15O9RO2'; 
fCHEXO2(i)=fCHEXO2(i)-1; fAPC10H15O9RO2(i)=fAPC10H15O9RO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1; 

%% APC10H15O10RO2 rxns
i=i+1;
Rnames{i} = 'APINAO2 + APC10H15O10RO2 = C20H32O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APINAO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fAPINAO2(i)=fAPINAO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H32O11(i)=fC20H32O11(i)+1; 

i=i+1;
Rnames{i} = 'APINBO2 + APC10H15O10RO2 = C20H32O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APINBO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fAPINBO2(i)=fAPINBO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H32O11(i)=fC20H32O11(i)+1; 

i=i+1;
Rnames{i} = 'APINCO2 + APC10H15O10RO2 = C20H32O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'APINCO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fAPINCO2(i)=fAPINCO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H32O11(i)=fC20H32O11(i)+1; 

i=i+1;
Rnames{i} = 'C107O2 + APC10H15O10RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C107O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC107O2(i)=fC107O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1; 

i=i+1;
Rnames{i} = 'C109O2 + APC10H15O10RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C109O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC109O2(i)=fC109O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1; 

i=i+1;
Rnames{i} = 'C106O2 + APC10H15O10RO2 = C20H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C106O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC106O2(i)=fC106O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O13(i)=fC20H30O13(i)+1; 

i=i+1;
Rnames{i} = 'C920CO3 + APC10H15O10RO2 = C20H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C920CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC920CO3(i)=fC920CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O13(i)=fC20H30O13(i)+1; 

i=i+1;
Rnames{i} = 'C108O2 + APC10H15O10RO2 = C20H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C108O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC108O2(i)=fC108O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O13(i)=fC20H30O13(i)+1; 

i=i+1;
Rnames{i} = 'PINALO2 + APC10H15O10RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'PINALO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fPINALO2(i)=fPINALO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1; 

i=i+1;
Rnames{i} = 'C96CO3 + APC10H15O10RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C96CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC96CO3(i)=fC96CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1; 

i=i+1;
Rnames{i} = 'C923CO3 + APC10H15O10RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C923CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC923CO3(i)=fC923CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1; 

i=i+1;
Rnames{i} = 'LIMAO2 + APC10H15O10RO2 = C20H32O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMAO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fLIMAO2(i)=fLIMAO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H32O11(i)=fC20H32O11(i)+1; 

i=i+1;
Rnames{i} = 'LIMALBO2 + APC10H15O10RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMALBO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fLIMALBO2(i)=fLIMALBO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1; 

i=i+1;
Rnames{i} = 'LIMCO2 + APC10H15O10RO2 = C20H32O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMCO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fLIMCO2(i)=fLIMCO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H32O11(i)=fC20H32O11(i)+1; 

i=i+1;
Rnames{i} = 'LIMALO2 + APC10H15O10RO2 = C20H32O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'LIMALO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fLIMALO2(i)=fLIMALO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H32O13(i)=fC20H32O13(i)+1; 

i=i+1;
Rnames{i} = 'LIMBO2 + APC10H15O10RO2 = C20H32O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMBO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fLIMBO2(i)=fLIMBO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H32O11(i)=fC20H32O11(i)+1; 

i=i+1;
Rnames{i} = 'LIMALAO2 + APC10H15O10RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LIMALAO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fLIMALAO2(i)=fLIMALAO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1; 

i=i+1;
Rnames{i} = 'NAPINAO2 + APC10H15O10RO2 = C20H31O10NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NAPINAO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNAPINAO2(i)=fNAPINAO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H31O10NO3(i)=fC20H31O10NO3(i)+1; 

i=i+1;
Rnames{i} = 'NAPINBO2 + APC10H15O10RO2 = C20H31O10NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NAPINBO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNAPINBO2(i)=fNAPINBO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H31O10NO3(i)=fC20H31O10NO3(i)+1; 

i=i+1;
Rnames{i} = 'BPINAO2 + APC10H15O10RO2 = C20H32O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'BPINAO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fBPINAO2(i)=fBPINAO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H32O11(i)=fC20H32O11(i)+1; 

i=i+1;
Rnames{i} = 'BPINBO2 + APC10H15O10RO2 = C20H32O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'BPINBO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fBPINBO2(i)=fBPINBO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H32O11(i)=fC20H32O11(i)+1; 

i=i+1;
Rnames{i} = 'BPINCO2 + APC10H15O10RO2 = C20H32O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'BPINCO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fBPINCO2(i)=fBPINCO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H32O11(i)=fC20H32O11(i)+1; 

i=i+1;
Rnames{i} = 'C918CO3 + APC10H15O10RO2 = C20H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C918CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC918CO3(i)=fC918CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H30O12(i)=fC20H30O12(i)+1; 

i=i+1;
Rnames{i} = 'NC102O2 + APC10H15O10RO2 = C20H29O12NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NC102O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNC102O2(i)=fNC102O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H29O12NO3(i)=fC20H29O12NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC101O2 + APC10H15O10RO2 = C20H29O11NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NC101O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNC101O2(i)=fNC101O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H29O11NO3(i)=fC20H29O11NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMO2 + APC10H15O10RO2 = C20H31O10NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NLIMO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNLIMO2(i)=fNLIMO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H31O10NO3(i)=fC20H31O10NO3(i)+1; 

i=i+1;
Rnames{i} = 'NLIMALO2 + APC10H15O10RO2 = C20H31O12NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NLIMALO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNLIMALO2(i)=fNLIMALO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H31O12NO3(i)=fC20H31O12NO3(i)+1; 

i=i+1;
Rnames{i} = 'NC91CO3 + APC10H15O10RO2 = C20H29O11NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NC91CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNC91CO3(i)=fNC91CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H29O11NO3(i)=fC20H29O11NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINAO2 + APC10H15O10RO2 = C20H31O10NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NBPINAO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNBPINAO2(i)=fNBPINAO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H31O10NO3(i)=fC20H31O10NO3(i)+1; 

i=i+1;
Rnames{i} = 'NBPINBO2 + APC10H15O10RO2 = C20H31O10NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NBPINBO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNBPINBO2(i)=fNBPINBO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC20H31O10NO3(i)=fC20H31O10NO3(i)+1; 

i=i+1;
Rnames{i} = 'C96O2 + APC10H15O10RO2 = C19H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C96O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC96O2(i)=fC96O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H30O11(i)=fC19H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C89CO3 + APC10H15O10RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C89CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC89CO3(i)=fC89CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C920O2 + APC10H15O10RO2 = C19H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C920O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC920O2(i)=fC920O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H30O12(i)=fC19H30O12(i)+1; 

i=i+1;
Rnames{i} = 'C97O2 + APC10H15O10RO2 = C19H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C97O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC97O2(i)=fC97O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H30O12(i)=fC19H30O12(i)+1; 

i=i+1;
Rnames{i} = 'C85CO3 + APC10H15O10RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C85CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC85CO3(i)=fC85CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C811CO3 + APC10H15O10RO2 = C19H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C811CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC811CO3(i)=fC811CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O13(i)=fC19H28O13(i)+1; 

i=i+1;
Rnames{i} = 'C921O2 + APC10H15O10RO2 = C19H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C921O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC921O2(i)=fC921O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H30O13(i)=fC19H30O13(i)+1; 

i=i+1;
Rnames{i} = 'C98O2 + APC10H15O10RO2 = C19H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C98O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC98O2(i)=fC98O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H30O13(i)=fC19H30O13(i)+1; 

i=i+1;
Rnames{i} = 'C922O2 + APC10H15O10RO2 = C19H30O14';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C922O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC922O2(i)=fC922O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H30O14(i)=fC19H30O14(i)+1; 

i=i+1;
Rnames{i} = 'C923O2 + APC10H15O10RO2 = C19H30O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C923O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC923O2(i)=fC923O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H30O11(i)=fC19H30O11(i)+1; 

i=i+1;
Rnames{i} = 'C924O2 + APC10H15O10RO2 = C19H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C924O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC924O2(i)=fC924O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H30O12(i)=fC19H30O12(i)+1; 

i=i+1;
Rnames{i} = 'C816CO3 + APC10H15O10RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C816CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC816CO3(i)=fC816CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'NORLIMO2 + APC10H15O10RO2 = C19H30O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NORLIMO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNORLIMO2(i)=fNORLIMO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H30O13(i)=fC19H30O13(i)+1; 

i=i+1;
Rnames{i} = 'LMKAO2 + APC10H15O10RO2 = C19H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LMKAO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fLMKAO2(i)=fLMKAO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H30O12(i)=fC19H30O12(i)+1; 

i=i+1;
Rnames{i} = 'LMKBO2 + APC10H15O10RO2 = C19H30O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'LMKBO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fLMKBO2(i)=fLMKBO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H30O12(i)=fC19H30O12(i)+1; 

i=i+1;
Rnames{i} = 'C926O2 + APC10H15O10RO2 = C19H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C926O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC926O2(i)=fC926O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O13(i)=fC19H28O13(i)+1; 

i=i+1;
Rnames{i} = 'C817CO3 + APC10H15O10RO2 = C19H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C817CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC817CO3(i)=fC817CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O13(i)=fC19H28O13(i)+1; 

i=i+1;
Rnames{i} = 'LMLKAO2 + APC10H15O10RO2 = C19H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'LMLKAO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fLMLKAO2(i)=fLMLKAO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O13(i)=fC19H28O13(i)+1; 

i=i+1;
Rnames{i} = 'LMLKBO2 + APC10H15O10RO2 = C19H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'LMLKBO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fLMLKBO2(i)=fLMLKBO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O13(i)=fC19H28O13(i)+1; 

i=i+1;
Rnames{i} = 'C823CO3 + APC10H15O10RO2 = C19H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C823CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC823CO3(i)=fC823CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O13(i)=fC19H28O13(i)+1; 

i=i+1;
Rnames{i} = 'C925O2 + APC10H15O10RO2 = C19H30O14';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C925O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC925O2(i)=fC925O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H30O14(i)=fC19H30O14(i)+1; 

i=i+1;
Rnames{i} = 'NOPINAO2 + APC10H15O10RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NOPINAO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNOPINAO2(i)=fNOPINAO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'NOPINBO2 + APC10H15O10RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NOPINBO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNOPINBO2(i)=fNOPINBO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'NOPINCO2 + APC10H15O10RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NOPINCO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNOPINCO2(i)=fNOPINCO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'NOPINDO2 + APC10H15O10RO2 = C19H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'NOPINDO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNOPINDO2(i)=fNOPINDO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O11(i)=fC19H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C918O2 + APC10H15O10RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C918O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC918O2(i)=fC918O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C9DCO2 + APC10H15O10RO2 = C19H26O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C9DCO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC9DCO2(i)=fC9DCO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H26O12(i)=fC19H26O12(i)+1; 

i=i+1;
Rnames{i} = 'C915O2 + APC10H15O10RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C915O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC915O2(i)=fC915O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C917O2 + APC10H15O10RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C917O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC917O2(i)=fC917O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C919O2 + APC10H15O10RO2 = C19H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C919O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC919O2(i)=fC919O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O13(i)=fC19H28O13(i)+1; 

i=i+1;
Rnames{i} = 'C914O2 + APC10H15O10RO2 = C19H26O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C914O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC914O2(i)=fC914O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H26O13(i)=fC19H26O13(i)+1; 

i=i+1;
Rnames{i} = 'C916O2 + APC10H15O10RO2 = C19H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C916O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC916O2(i)=fC916O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O13(i)=fC19H28O13(i)+1; 

i=i+1;
Rnames{i} = 'C88CO3 + APC10H15O10RO2 = C19H26O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C88CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC88CO3(i)=fC88CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H26O13(i)=fC19H26O13(i)+1; 

i=i+1;
Rnames{i} = 'C87CO3 + APC10H15O10RO2 = C19H26O14';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C87CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC87CO3(i)=fC87CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H26O14(i)=fC19H26O14(i)+1; 

i=i+1;
Rnames{i} = 'C822CO3 + APC10H15O10RO2 = C19H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C822CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC822CO3(i)=fC822CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H28O12(i)=fC19H28O12(i)+1; 

i=i+1;
Rnames{i} = 'NLMKAO2 + APC10H15O10RO2 = C19H29O11NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NLMKAO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNLMKAO2(i)=fNLMKAO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC19H29O11NO3(i)=fC19H29O11NO3(i)+1; 

i=i+1;
Rnames{i} = 'C85O2 + APC10H15O10RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C85O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC85O2(i)=fC85O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C89O2 + APC10H15O10RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C89O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC89O2(i)=fC89O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C86O2 + APC10H15O10RO2 = C18H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C86O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC86O2(i)=fC86O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O12(i)=fC18H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C811O2 + APC10H15O10RO2 = C18H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C811O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC811O2(i)=fC811O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O12(i)=fC18H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C810O2 + APC10H15O10RO2 = C18H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C810O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC810O2(i)=fC810O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O12(i)=fC18H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C812O2 + APC10H15O10RO2 = C18H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C812O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC812O2(i)=fC812O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O13(i)=fC18H28O13(i)+1; 

i=i+1;
Rnames{i} = 'C813O2 + APC10H15O10RO2 = C18H28O14';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C813O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC813O2(i)=fC813O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O14(i)=fC18H28O14(i)+1; 

i=i+1;
Rnames{i} = 'C729CO3 + APC10H15O10RO2 = C18H26O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C729CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC729CO3(i)=fC729CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H26O12(i)=fC18H26O12(i)+1; 

i=i+1;
Rnames{i} = 'C816O2 + APC10H15O10RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C816O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC816O2(i)=fC816O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C817O2 + APC10H15O10RO2 = C18H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C817O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC817O2(i)=fC817O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O12(i)=fC18H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C826O2 + APC10H15O10RO2 = C18H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C826O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC826O2(i)=fC826O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O13(i)=fC18H28O13(i)+1; 

i=i+1;
Rnames{i} = 'C822O2 + APC10H15O10RO2 = C18H28O11';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C822O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC822O2(i)=fC822O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O11(i)=fC18H28O11(i)+1; 

i=i+1;
Rnames{i} = 'C818O2 + APC10H15O10RO2 = C18H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C818O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC818O2(i)=fC818O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O13(i)=fC18H28O13(i)+1; 

i=i+1;
Rnames{i} = 'C823O2 + APC10H15O10RO2 = C18H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C823O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC823O2(i)=fC823O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O12(i)=fC18H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C819O2 + APC10H15O10RO2 = C18H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C819O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC819O2(i)=fC819O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O13(i)=fC18H28O13(i)+1; 

i=i+1;
Rnames{i} = 'C727CO3 + APC10H15O10RO2 = C18H26O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C727CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC727CO3(i)=fC727CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H26O13(i)=fC18H26O13(i)+1; 

i=i+1;
Rnames{i} = 'C731CO3 + APC10H15O10RO2 = C18H26O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C731CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC731CO3(i)=fC731CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H26O13(i)=fC18H26O13(i)+1; 

i=i+1;
Rnames{i} = 'C824O2 + APC10H15O10RO2 = C18H28O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C824O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC824O2(i)=fC824O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O12(i)=fC18H28O12(i)+1; 

i=i+1;
Rnames{i} = 'C820O2 + APC10H15O10RO2 = C18H26O14';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C820O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC820O2(i)=fC820O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H26O14(i)=fC18H26O14(i)+1; 

i=i+1;
Rnames{i} = 'C825O2 + APC10H15O10RO2 = C18H28O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C825O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC825O2(i)=fC825O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O13(i)=fC18H28O13(i)+1; 

i=i+1;
Rnames{i} = 'C732CO3 + APC10H15O10RO2 = C18H26O14';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C732CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC732CO3(i)=fC732CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H26O14(i)=fC18H26O14(i)+1; 

i=i+1;
Rnames{i} = 'C8BCO2 + APC10H15O10RO2 = C18H28O10';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C8BCO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC8BCO2(i)=fC8BCO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H28O10(i)=fC18H28O10(i)+1; 

i=i+1;
Rnames{i} = 'C88O2 + APC10H15O10RO2 = C18H26O12';
k(:,i) = kROORc.*fROOR;
Gstr{i,1} = 'C88O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC88O2(i)=fC88O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H26O12(i)=fC18H26O12(i)+1; 

i=i+1;
Rnames{i} = 'C718CO3 + APC10H15O10RO2 = C18H26O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C718CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC718CO3(i)=fC718CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H26O13(i)=fC18H26O13(i)+1; 

i=i+1;
Rnames{i} = 'C87O2 + APC10H15O10RO2 = C18H26O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C87O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC87O2(i)=fC87O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H26O13(i)=fC18H26O13(i)+1; 

i=i+1;
Rnames{i} = 'C721CO3 + APC10H15O10RO2 = C18H26O13';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'C721CO3'; Gstr{i,2} = 'APC10H15O10RO2'; 
fC721CO3(i)=fC721CO3(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H26O13(i)=fC18H26O13(i)+1; 

i=i+1;
Rnames{i} = 'NC826O2 + APC10H15O10RO2 = C18H27O12NO3';
k(:,i) = kROORd.*fROOR;
Gstr{i,1} = 'NC826O2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fNC826O2(i)=fNC826O2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fC18H27O12NO3(i)=fC18H27O12NO3(i)+1; 

i=i+1;
Rnames{i} = 'BUT2OLO2 + APC10H15O10RO2 = BUT2OLdimer';
k(:,i) = kROORc.*fROOR1;
Gstr{i,1} = 'BUT2OLO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fBUT2OLO2(i)=fBUT2OLO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fBUT2OLdimer(i)=fBUT2OLdimer(i)+1; 

i=i+1;
Rnames{i} = 'CHEXO2 + APC10H15O10RO2 = CHEXdimer';
k(:,i) = kROORc.*fROOR2;
Gstr{i,1} = 'CHEXO2'; Gstr{i,2} = 'APC10H15O10RO2'; 
fCHEXO2(i)=fCHEXO2(i)-1; fAPC10H15O10RO2(i)=fAPC10H15O10RO2(i)-1; fCHEXdimer(i)=fCHEXdimer(i)+1; 

%END OF REACTION LIST