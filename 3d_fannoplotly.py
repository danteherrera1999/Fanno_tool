import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
import plotly.graph_objects as go
#Constants
in2m = 0.0254
psi2pa = 6894.75729
m2ft = 1 / (12 * in2m)
gamma_0 = 1.299
gamma = gamma_0
R_bar = 8314.46261815324
gpm2cfm = 0.13368055363
#BASE EQUATIONS
Re = lambda rho, V, D, mu: rho * V * D / mu
dff = lambda x,Re,epsilon,D: 1/np.sqrt(x[0]) + 2 * np.log10((epsilon / (D *  3.7)) + (2.51 / (Re * np.sqrt(x[0]))))
flod = lambda M,gamma=gamma_0,n=0: (gamma + 1)/ 2 / gamma * np.log(((gamma + 1)/2)/(1 + (gamma - 1) / 2 * M ** 2)) - (1 - M ** -2) / gamma - ((gamma + 1) / 2 / gamma) * np.log(M ** -2) - n
tots = lambda M,gamma=gamma_0: (gamma + 1) / (2 + (gamma - 1) * M ** 2)
pops = lambda M,gamma=gamma_0: ((gamma + 1)/(2 + (gamma - 1) * M ** 2)) ** (1/2) / M
rors = lambda M,gamma=gamma_0: ((2 + (gamma - 1) * M ** 2) / (gamma + 1)) ** (1/2) / M
pnopns = lambda M,gamma=gamma_0 : (((2 + (gamma - 1) * M ** 2) / (gamma + 1) ) ** ((gamma + 1) / (2 * (gamma - 1) )) ) / M
Mfflod = lambda flod_r,gamma=gamma_0: fsolve(flod,.01,args=(gamma,flod_r))[0]
MfLm = lambda L,D,f: Mfflod(f * L / D)
mdot2scfm = lambda m_dot, rho_STP: m_dot*2.205*60/(0.062428*rho_STP)
scfm2mdot = lambda SCFM: rho_STP * SCFM * ft3pm2m3ps
#STP CONDITIONS
WF = "CH4"
T_STP = 293.15
P_STP = 14.7 * psi2pa
rho_air_STP = PropsSI("D","P",P_STP,"T",T_STP,"Air")
rho_STP = PropsSI("D","P",P_STP,"T",T_STP,WF)
MW = 16.043
MW_air = 28.9647
epsilon = .045 / 1000
#DEFINE EXIT CONDITIONS
T_e = T_STP
P_e = P_STP
rho_e = PropsSI("D","P",P_e,"T",T_e,WF)
mu_e = PropsSI("V","P",P_e,"T",T_e,WF)
#ITERATION PARAMETERS
f_rtgs = .05 #INITIAL GUESS OF FRICTION FACTOR
n = 3 #NUMBER OF ITERATIONS
rf = 10 ** 7 #ROUND FACTOR
c_p = .878
P_upstream = 425
g = 40 #GRID PARAMETER (N>2 USE WISELY)
#DOES IT CONVERGE ?!?(yes :3)
def fiter(f,L,D,M_e,T_s,P_s):
    L_max2 = flod(M_e) * D / f
    M_inlet = MfLm(L + L_max2,D,f)
    #INLET CONDITIONS
    T_inlet = tots(M_inlet) * T_s
    P_inlet = pops(M_inlet) * P_s
    mu_inlet = PropsSI("V","P",P_inlet,"T",T_inlet,WF)
    rho_inlet = PropsSI("D","P",P_inlet,"T",T_inlet,WF)
    a_inlet = np.sqrt(R_bar / MW * gamma * T_inlet)
    V_inlet = a_inlet * M_inlet
    Re_inlet = Re(rho_inlet,V_inlet,D,mu_inlet)
    A = np.pi * D ** 2 / 4
    Q = A * V_inlet
    SCFM = Q * P_inlet / P_s * T_s / T_inlet * np.sqrt(MW/MW_air)
    return M_inlet, fsolve(dff,.01,args=(Re_inlet,epsilon,D)),SCFM
f = f_rtgs
#MAIN ITERATION LOOP FOR 1 EXIT MACH NUMBER
def run(f,L,D,M_e):
    P_s = P_e / pops(M_e)
    T_s = T_e / tots(M_e)
    for D in Ds:
        for L in Ls:
            for i in range(n):
                M, f,SCFM = fiter(f,L,D,M_e,T_s,P_s)
            Ps.append(pops(M)*P_s/psi2pa)
            SCFMS.append(SCFM)
Ds = np.linspace(1,11,g) * in2m
Ls = np.linspace(0,1000/m2ft,g)
M_es = np.arange(.5,1.05,.1)
#PLOTTING STUFF
X, Y = np.meshgrid(Ls*m2ft,Ds/in2m)
colorscale = [[0,"blue"],[300,"red"]]
surfaces = []
for M_e in M_es:
    Ps = []
    SCFMS = []
    run(f_rtgs,Ls,Ds,M_e)
    Ps = np.reshape(np.array(Ps),(g,g))
    SCFMS = np.reshape(np.array(SCFMS), (g, g))
    surfaces.append(go.Surface(x=X,y=Y,z=Ps,cmax=300,cmin=0,colorbar={'title':"Pressure (PSI)"},colorscale='plasma',name=f"M_e: {M_e}"))
layout = go.Layout(title='Effect of Pipe Diameter and Length on Inlet Pressure',title_x=.5,
                legend=dict(
                orientation="h")
                )
fig = go.Figure(data=surfaces,layout=layout)
fig.update_scenes(xaxis_title_text='Pipe Length (ft)',
                  yaxis_title_text='Pipe Diameter (In)',
                  zaxis_title_text='Inlet Pressure (PSI)')

fig.show()
