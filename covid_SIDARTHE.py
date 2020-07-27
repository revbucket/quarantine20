""" Data from this April Nature paper: 
https://www.nature.com/articles/s41591-020-0883-7#data-availability


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation horizon: CAN BE MODIFIED AT ONE'S WILL PROVIDED THAT IT IS AT
% LEAST EQUAL TO THE NUMBER OF DAYS FOR WHICH DATA ARE AVAILABLE
Orizzonte = 350;

% Plot yes/no: SET TO 1 IF PDF FIGURES MUST BE GENERATED, 0 OTHERWISE
plotPDF = 0;

% Time-step for Euler discretisation of the continuous-time system
step=0.01;
"""

# Transmission rate due to contacts with UNDETECTED asymptomatic infected
alfa=0.57
# Transmission rate due to contacts with DETECTED asymptomatic infected
beta=0.0114
# Transmission rate due to contacts with UNDETECTED symptomatic infected
gamma=0.456
# Transmission rate due to contacts with DETECTED symptomatic infected
delta=0.0114

# Detection rate for ASYMPTOMATIC
epsilon=0.171
# Detection rate for SYMPTOMATIC
theta=0.3705

# Worsening rate: UNDETECTED asymptomatic infected becomes symptomatic
zeta=0.1254
# Worsening rate: DETECTED asymptomatic infected becomes symptomatic
eta=0.1254

# Worsening rate: UNDETECTED symptomatic infected develop life-threatening
# symptoms
mu=0.0171
# Worsening rate: DETECTED symptomatic infected develop life-threatening
# symptoms
nu=0.0274

# Mortality rate for infected with life-threatening symptoms
tau=0.01

# Recovery rate for undetected asymptomatic infected
lambda_ =0.0342
# Recovery rate for detected asymptomatic infected
rho=0.0342
# Recovery rate for undetected symptomatic infected
kappa=0.0171
# Recovery rate for detected symptomatic infected
xi=0.0171
# Recovery rate for life-threatened symptomatic infected
sigma=0.0171




"""
SIDHARTE MODEL != SIR 
(but we can model it kinda similarly with...)
S-> I rate := alpha + beta + gamma + delta 
I-> R rate: lambda + rho + kappa + xi + sigma + tau 
"""

SIR_infect = alfa + beta + gamma +delta 
SIR_recover = lambda_ + rho + kappa + xi + sigma + tau 

print("SIR INFECT", SIR_infect)
print("SIR RECOVER", SIR_recover)