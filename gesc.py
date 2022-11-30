## Direct translation of Marine's GESC code from MATLAB to Python

## AIR  - SOLID with 2 solid layers over half space
# Example of a simple gradient over half space. 
# Marine Denolle (09/2014)

import numpy as np

# Create time/frequency structure FT
def make_ft(win,rate,typ):
# input
# choose between time (typ='t') or frequency ('f'):
#     enter win: window length (in s or in Hz)
#     enter rate : sampling rate (in 1/s or 1/Hz)
#  choose either the frequencies to solve the eigen problem at,
#  AND/OR the time domain to build seismograms
#  Marine Denolle (04/2014)

    if typ=='t': # if TIME DOMAIN
        Fnyq = 1/(2*rate)
        nwin = np.floor(win/rate)
        omega = 2*np.pi*np.linspace(0,Fnyq,int(nwin/2+1))
        FT = {'twin': win,        # (in s) time series length
              'dt': rate,         # (in 1/s) sampling rate
              'Nsps': 1/rate, # number of sample per seconds 
              'nwin': nwin, # # points in time series
              'Fnyq': Fnyq, # Nyquist frequency
              'omega': omega, # angular frequency vector
              'df': omega[1]/(2*np.pi),
              'tt': np.arange(0,win-rate,rate),
              'Fmax':Fnyq}
#         FT.tt=(0:FT.dt:FT.twin-FT.dt)
    elif typ=='f':
        print('not implemented yet')
#         FT.Fnyq = win # (in Hz) Fmax
#         FT.df = rate # (in 1/Hz) sampling rate
#         FT.dt = 1/(2*FT.Fnyq) # (in 1/s) sampling rate of hypothetical time series
#         FT.nwin = 2*floor(FT.Fnyq/FT.df)-1 # # of hypothetical time series
#         FT.omega=2*pi*linspace(0,FT.Fnyq,floor(FT.Fnyq/FT.df)) # angular frequency vector
#         FT.Nsps = 1/FT.dt
#     FT['Fmax']=FT['Fnyq'] # by default, the max frequency to look at is the Nyquist frequency

    return FT

## Create Medium structure MED using layers
def make_layers(H,Dmax,alpha,beta,rho,typ):
    Nl=len(H)
    ## input
    # H: (in km) vector of interface depths starting with depth = 0 km (at the
    # surface) and positive downward
    # Dmax: input thickness of half space: bottom layer may be taken thick to
    # simulate a half space: take it large for broadband caluclations
    # alpha: (in km/s) vector of P-wave speed (Vp) of length *length(H)-1*
    # beta: (in km/s) vector of S-wave speed (Vs) of length *length(H)-1*
    # rho: (in kg/dm3) vector of density of length *length(H)-1*
    # typ: vector of {'cst'} (homogeneous layers), {'lin'} (constant gradient layers)

    ## output
    # MED: structure with medium properties. GESC requires organizing medium
    # from bottom to top, and for elastic properties *within* layers to match
    # properly boundary conditions

    # the GESC code asks for the elastic properties from the bottom of the
    # MEDium to the top:
    # MED is a vector (length = Nl number of layers) of structures
    # MED[0] contains depth, Vs, Vp, Rho and depth of interfaces and solid-water
    # flag for entire MEDium
    # MED[0] corresponds to the half space
    # MED(2:Nl) contains depth, Vs, Vp, Rho *within* the layers, from bottom to
    # top

    # Note: to choose a Poisson medium, poiss=1 and alpha = sqrt(3)beta
    #       to choose homogeneous layers, use typ = 'cst'
    #       to choose constant gradient layers, use typ = 'lin'
    # Marine Denolle (04/2014)


    # check if medium is Poisson (i.e. if Vp = sqrt(3)*Vs)
    MED = []
    MED.append({})
    MED[0]['poiss']=0 # not a Poisson medium
    
    if all(alpha==np.sqrt(3)*beta):  MED[1]['poiss'] = 1

    # check if all layers are solid (i.e. if Vs != 0)
    MED[0]['solid']=np.zeros(Nl)
    MED[0]['solid'][0]=1

    for ii in range(len(beta)):
        if beta[ii]!=0: 
            MED[0]['solid'][Nl-ii-1]=1


    # create medium in vectors
    z=np.array([])
    a=np.array([])
    b=np.array([])
    r=np.array([])
    alpha = np.hstack([alpha, alpha[-1] ])
    beta = np.hstack([beta,   beta[-1]  ])
    rho = np.hstack([rho, rho[-1]])
    for i in range(1,Nl):
        if typ[i-1]=='cst':
            z = np.hstack([np.linspace(H[i],H[i-1],50), z]) # depth vector
            a = np.hstack([alpha[i-1]*np.ones(50), a])  # alpha / P wave velocity vector
            b = np.hstack([beta[i-1]*np.ones(50), b])   # beta / S wave velocity vector
            r = np.hstack([rho[i-1]*np.ones(50), r])   # rtho /  density vector
        elif typ[i-1]=='lin':
            z = np.hstack([np.linspace(H[i],H[i-1],50), z])         # depth vector
            a = np.hstack([np.linspace(alpha[i],alpha[i-1],50), a]) # alpha / P wave velocity vector
            b = np.hstack([np.linspace(beta[i],beta[i-1],50), b])   # beta / S wave velocity vector
            r = np.hstack([np.linspace(rho[i],rho[i-1],50), r])     # rho /  density vector
        else:
            print('Enter type of layer poperly')
            return []

    # remove overlapping point (2 points have the same depth H1, H1)
    trash,ib = np.unique(z, return_index=True)
     
#     # outpout data in MED
    MED[0]['z'] = z[ib]
    MED[0]['alpha'] = a[ib]
    MED[0]['beta'] = b[ib]
    MED[0]['rho'] = r[ib]
    MED[0]['inter'] = H
    MED[0]['Nl'] = Nl 
               
    # break into layers: elastic properties have to be contained within layer
    temp_layers = []
    for i in range(1,MED[0]['Nl']):
        
        ik= np.where(      (MED[0]['inter'][i-1] < MED[0]['z'])\
                         & (MED[0]['inter'][i]   > MED[0]['z']) )[0]
        
        new_layer = {
            'betal':  MED[0]['beta'][ik], #[::-1],
            'alphal': MED[0]['alpha'][ik], #[::-1], 
            'rhol' :  MED[0]['rho'][ik], #[::-1],
            'zz' :    MED[0]['z'][ik], #[::-1],
            'typ' : typ[i-1]
        }
        
        
        new_layer['zz'][0]  = np.min( [MED[0]['inter'][i], MED[0]['inter'][i-1]] )
        new_layer['zz'][-1] = np.max( [MED[0]['inter'][i], MED[0]['inter'][i-1]] )
        
        temp_layers.append(new_layer)
        
    for layer in reversed(temp_layers):
        MED.append(layer)

    MED[0]['inter'] = MED[0]['inter'][::-1]
    MED[0]['Dmax'] = Dmax

    return MED

from numpy.matlib import repmat

def cheb(N):
    # CHEB  compute D = differentiation matrix, x = Chebyshev grid
    if N==0: return 0,1
    x = np.cos(np.pi*(np.arange(N+1))/N)
    x = np.array([x]).T
    c = np.hstack([2, np.ones(N-1), 2]) * np.power(-1,np.arange(N+1)) #
    c = np.array([c])
    X = repmat(x,1,N+1);
    dX = X-X.T;                  
    D  = np.matmul(c.T,(1/c)) / (dX+(np.eye(N+1)));      # off-diagonal entries
    D  = D - np.diag(sum(D.T));               # diagonal entries
    return D,x

# Subroutine to solve eigenpropblem given a medium and frequency
# Marine Denolle (04/10/14)
 
from scipy.interpolate import interp1d
    
def gesc(MED_GESC,FT,NR,max_mode,freq):

    ## OUTPUTS
    # bigC: structure that contains the medium properties per layer, the
    # eigenfunctions, the dispersion curves etc etc

    ## INPUTS
    # MED:medium structure
    # FT: time/frequency structure
    # NR: structure of variation of number of collocation points with
    # wavelength (contains for upper layer
    # NR : number of points scheme: see Section 2.3 on resolution 
    #     NR['typ'] = 'cst','lin' or 'exp' for type of variation of #points with
    #     thickness
    #     NR.rate = 'cst' (constant slope), 'rec' (recursive and depends on
    #     previous variations of the eigen functions)
    #     NR.lambd = depths at which change scheme to lower resolution number 
    #     NR.amp(1:length(NR.lambd)) = slope value (see Figure XX from manual)
    #     NR.N1(1:length(NR.lambd)) = lower bound # points
    #     NR.N2(1:length(NR.lambd)) = upper bound # points

    # max_mode: maximum number of modes to look at
    # freq(optional): user may ask to output eigenfunctions at frequencies stored in
    # the freq vector.
    # bigC(1:min(length(FT.omega),floor(FT.Fmax/FT.df+1)))=struct('a',[]);

    
#     bigC(1:min(length(FT.omega),floor(FT.Fmax/min(FT.df)+1)))=struct('a',[]);
    C = [{}]
    # calculate for each frequency, up to Nyquist or Fmax.
    if np. floor(FT['Fmax']/FT['df'])>len(FT['omega']):
        # assumes len(FT['df'])==1 
        error('!! Choose Fmax < Fnyq !!')
    nerr=1

    nf = len(FT['omega'])
    for ifreq in range(nf):
        if FT['omega'][ifreq]==0: continue

        
        MED=MED_GESC

        # If most of the eigenfunction variations is trapped within the upper
        # layer, we split that layer to reduce the number of points required to
        # solve the strong variations in the layer. Dichotomy process!

        if (ifreq==0) | (nerr==1):
            C1 = [{}]
            C1[0]['cl'] = []
            C1[0]['cl'].append ( max(MED[0]['beta']) )
            C1[0]['cr'] = []
            C1[0]['cr'].append ( max(MED[0]['beta']) )
            C1[0]['Nl'] = MED[0]['Nl']
            HH =  max(MED_GESC[0]['inter'])
        else:
            ikk=find(C1[0]['uz'][:,1]<=1E-2*abs(C1[0]['uz'][end,1]))
            HH=C1[0]['zz'](ikk(end))

            
        if (ifreq>0) & (HH< MED_GESC[0]['inter'][-2]):
            C.append({})

            if MED_GESC[0]['typ'] == 'cst':  # if this a layered medium with homoheneous layers
                for i in range(len(MED_GESC[0]['inter'])):
                    ikk=np.argmin(abs(MED_GESC[0]['z']-MED_GESC[0]['inter'][i]))
                    aa[i] = MED_GESC[0]['alpha'][ikk]
                    bb[i] = MED_GESC[0]['beta'][ikk]
                    rr[i] = MED_GESC[0]['rho'][ikk]

                aa[i+1] = aa[i]
                bb[i+1] = bb[i]
                rr[i+1] = rr[i]
                
                aa=fliplr(aa);
                bb=fliplr(bb);
                rr=fliplr(rr);
                
                MED[0]['inter'] = [0, HH, fliplr(MED_GESC[0]['inter'][1:-2])]
                MED= make_layers(MED[0]['inter'],MED_GESC[0]['Dmax'],aa,bb,rr,MED_GESC[0]['typ']);

            elif MED_GESC[0]['typ'] == 'lin':  # if you have gradients
                i=MED_GESC[0]['Nl']
                MED[0]['inter']=[MED_GESC[0]['inter'][1:-2], HH, 0];
                MED[0]['Nl']=len(MED[0]['inter']);
                for i in range(MED[0]['Nl']-1,MED[0]['Nl']):
                    MED[i]['zz']     = linspace(MED[0]['inter'](i-1),MED[0]['inter'][i],30)
                    
                    MED[i]['alphal'] = interp1d(MED_GESC[0]['z'],
                                               MED_GESC[0]['alpha'],
                                               MED[i]['zz'],'linear')
                    
                    MED[i]['betal']  = interp1d(MED_GESC[0]['z'],
                                               MED_GESC[0]['beta'], 
                                               MED[i]['zz'],'linear')
                    
                    MED[i]['rhol']   = interp1d(MED_GESC[0]['z'],
                                               MED_GESC[0]['rho'],  
                                               MED[i]['zz'],'linear')
                    
            print([f'Added a layer because first wavelength is {HH}'  \
                   f"km compared to {MED[0]['inter'](end-1)} km"])

            
        C[0]['Nl'] = MED[0]['Nl']
        C[0]['omega']=FT['omega'][ifreq]
        C[0]['Ntot'] = 0
        C[0]['solid']=MED[0]['solid']

        T=2*np.pi/C[0]['omega'];
        print(f"Period: {T} s and Frequency {C[0]['omega']/(2*np.pi)} Hz")


        for i in range(C[0]['Nl']):
             # find appropriate number of points for each layer   
            if i==0: #halfspace
                C[i]['H'] = MED[0]['Dmax']-MED[0]['inter'][0] # thickness
            else:
                C.append({})
                C[i]['H'] = -MED[0]['inter'][i] + MED[0]['inter'][i-1]  # thickness


            if NR['typ']=='typ':
                LL = C[i]['H']/(C1[0]['cl'][0]*T) # wavelength
                junk,ii=np.where( MED[0]['inter'][i]/(C1[0]['cl'][0]*T) <= NR.lambd) #? 1e4
                junk,kk=np.min(np.abs(LL-NR['x']))
                C[i]['N'] = np.floor(NR[ii[0]]['npts'][kk])

            if NR['typ']=='rec':
                mgrd=1
                if ifreq>1:
                    # find peak gradient in that layer
                    if (C[0]['Nl']==C1[0]['Nl'] | i < C1[0]['Nl']) & nerr==0:
                        for imode in range(1): #?
                            n1=range(C1[i]['nn'][0]+1,C1[i]['nn'](2))
                            u1=C1[0]['ux'][n1,imode]
                            u2=C1[0]['uy'][n1,imode]
                            u3=C1[0]['uz'][n1,imode]
                            dr1 = max(abs(C1[i]['D']*u1));
                            dl1 = max(abs(C1[i]['D']*u2));
                            dr2 = max(abs(C1[i]['D']*u3));
                            grd[imode] = np.mean([dr1, dr2, dl1])
                        mgrd=mean(grd)
                C[i]['N'] = min(max(
                                np.floor( NR['N1']*mgrd*C[i]['H']),NR['N1']),NR['N2'])

            # now you have H and N, define Cheb D
            D, x = cheb(C[i]['N']-1)     # differentiation matrix
            C[i]['D'] = 2/(C[i]['H']) * D
            C[i]['z']= (C[i]['H']/2*(x+1)).T[0];        # all that normalized - deep layer


            # fill in the zz, alpha, beta, rho and derived elastic parameters for
            # each layer and at each depth:
            C[i]['alpha'] = np.zeros(C[i]['N'])
            C[i]['beta'] = np.zeros(C[i]['N'])
            C[i]['rho'] = np.zeros(C[i]['N'])
            C[i]['mu'] = np.zeros(C[i]['N'])
            C[i]['lambd'] = np.zeros(C[i]['N'])
            for k in range(C[i]['N']):

                if i==0:  # for half space "layer"
                    beta=max(MED[0]['beta'])
                    alpha=max(MED[0]['alpha'])
                    rho=max(MED[0]['rho'])
                else:
#                     junk,ii=np.unique(MED[i]['zz'],'stable'))
                    ii = np.sort(np.unique(MED[i]['zz'],return_index=True)[1])

                    b=MED[i]['betal'][ii]
                    a=MED[i]['alphal'][ii]
                    r=MED[i]['rhol'][ii]
                    zz = MED[0]['inter'][i]+C[i]['z'][k]

                    f = interp1d(MED[i]['zz'][ii],b)
                    beta = f(zz)
                    
                    f = interp1d(MED[i]['zz'][ii],r)
                    rho = f(zz)
                    
                    f = interp1d(MED[i]['zz'][ii],a)
                    alpha = f(zz)

                C[i]['alpha'][k]    = alpha
                C[i]['beta'][k]     = beta
                C[i]['rho'][k]      = rho
                C[i]['mu'][k]       = rho*beta**2
                C[i]['lambd'][k]   = rho*alpha**2 - 2*C[i]['mu'][k]


            C[i]['lambdmu']=C[i]['lambd']+2*C[i]['mu']
            
            if i==0:
                C[0]['beta1D'] = C[0]['beta']
                C[0]['alpha1D'] = C[0]['alpha']
                C[0]['rho1D'] = C[0]['rho']
                C[0]['zz'] = C[0]['z'] + MED[0]['inter'][0]
            else:
                C[0]['beta1D']  = np.hstack([C[0]['beta1D']  , C[i]['beta']  ])
                C[0]['alpha1D'] = np.hstack([C[0]['alpha1D'] , C[i]['alpha'] ])
                C[0]['rho1D']   = np.hstack([C[0]['rho1D']   , C[i]['rho']   ])
                C[0]['zz']      = np.hstack([C[0]['zz']      , C[i]['z'] + MED[0]['inter'][i]])

            C[0]['Ntot'] = C[0]['Ntot'] + C[i]['N']
#             Ntot[ifreq]=C[0]['Ntot'] # Never called?

        C[0]['mu1D'] = C[0]['rho1D']*C[0]['beta1D']**2
        C[0]['lambd1D'] = C[0]['rho1D']*C[0]['alpha1D']**2 - 2*C[0]['mu1D']
        # if poiss == 1;C[0].lambd1D = C[0].mu1D;C[0].beta1D=1/sqrt(3)*C[0].alpha1D;end
        C[0]['lambdmu1D'] = C[0]['lambd1D'] + 2*C[0]['mu1D']

        if ifreq>1 & C[0]['Ntot'] > 500:  return
        
        
        '''
        solve eigenproblem
        '''
        
        # Rayleigh waves:
        A = RayInnerMat(C,C[0]['omega'])
        A,B = RayBC(C,A)
        C,nerr = EigRW(C,A,B,max_mode)
        if (nerr==1): continue

        '''    
        # Love waves:
        A = LoveInnerMat(C,C[0].omega)
        A,B = LoveBC(C,A);
        C,nerr = EigLW(C,A,B,max_mode)
        if (nerr==1): continue
        C[0]['Nmode'] = min([max_mode, len(C[0].kl), len(C[0].kr)])
        '''

        # get group velocity and integrals:
        C=get_integrals_sw(C)

        C1=C
        bigC[ifreq] = struct('a',C)
        junk,ib=np.unique(C[0].zz);
        uyi[ifreq] = interp1interp1d(C[0].zz(ib),C[0].uy(ib),5,'linear');

    return bigC 

from numpy import diag
def RayInnerMat(C,omega,verbose=False):
    
    Nl = C[0]['Nl']   # number of layers
    A =  np.zeros((4*C[0]['Ntot'],4*C[0]['Ntot']))
    if verbose: print(f"A.shape={A.shape}")
    n = 0
    
    for i in range(Nl):
        if verbose: print(f"C[{i}]['N'] = {C[i]['N']}")
        I = slice(n, n + 4*C[i]['N'])
            
        zzz = np.zeros( (C[i]['N'],C[i]['N']) )
        a_01 = -diag(C[i]['lambd']/C[i]['lambdmu'])*C[i]['D']
        a_02 = diag(1/C[i]['lambdmu'])
        a_10 = C[i]['D']
        a_13 = -diag(1/C[i]['mu'])
        a_20 = diag(C[i]['rho'])*(omega*omega.T)
        a_31 = -diag(C[i]['rho'])*(omega*omega.T) - C[i]['D']*diag( 4*C[i]['mu']*\
                    (C[i]['lambd']+C[i]['mu'])/C[i]['lambdmu'])*C[i]['D']
        a_32 = -C[i]['D']*diag(C[i]['lambd']/C[i]['lambdmu'])
        
        if verbose: print(f"zzz.shape={zzz.shape}")
        test = np.vstack([ # bottom layer      
                    np.hstack([zzz  , a_01 , a_02 , zzz ]),
                    np.hstack([a_10 , zzz  , zzz  , a_13]),
                    np.hstack([a_20 , zzz  , zzz  , a_10]),
                    np.hstack([zzz  , a_31 , a_32 , zzz])
                ]) 
        if verbose: print(f"test.shape={test.shape}")
        if verbose: print(f"A[I,I].shape={A[I,I].shape}")
        A[I,I] = test
        n = n + 4*C[i]['N']
        
    return A

def RayBC(C,A):
    Nl = C[0]['Nl'];   # number of layers
    B  = np.eye(A.shape[0])

    '''
    Bottom BC's (remove two above)
    '''
    
    # Uz =0 bottom
    A[C[0]['N'],:] = 0
    A[C[0]['N'],1] = 1
    B[C[0]['N'],:] = 0

    # Ux = 0
    A[0,         :] = 0
    A[0, C[0]['N']] = 1
    B[0,         :] = 0

    # for each interface (replace one above and one below the interface)
    #disp('series 1')
    n = -1 # I think this single change fixes most mATLAB->python indexing
    for ii in range(Nl-1):
        Nii = C[ii]['N']
        n = n + 4*Nii
        j = n - Nii

        # Continuity conditions
        # stress sigma_zz
        
        A[j,:] = 0
        B[j,:] = 0
        A[j,n-3*Nii:n-2*Nii] = -C[ii  ]['lambdmu'][-1] * C[ii]['D'][-1, :]
        B[j,n-3*Nii]         = +C[ii  ]['lambd']  [-1]
        A[j,n+Nii:n+2*Nii]   = +C[ii+1]['lambdmu'][ 0] * C[ii]['D'][1  ,:]
        B[j,n+1]             = -C[ii+1]['lambd']  [ 0]

        # Continuity conditions on sigma_xz (R4)
        A[n,:] =  0
        A[n,n] = -1
        A[n,n+3*C[ii+1]['N']+1] = 1
        B[n,:] =  0


        # continuous r1 (ux)
        A[n+1,            :]    =  0
        A[n+1,n-3*C[ii  ]['N']] = -1
        A[n+1,n+1]              =  1
        B[n+1          ,  :]    =  0

        # Continuity conditions on uz (r2)
        A[n+1+C[ii+1]['N'],          :] =  0
        A[n+1+C[ii+1]['N'],n-2*C[ii]['N']] = -1
        A[n+1+C[ii+1]['N'],n+1+C[ii+1]['N']] =  1
        B[n+1+C[ii+1]['N'],          :] =  0


        # if the interface is water-solid, free shear condition
#         if max(C(ii+1).mu)==0 : # Never called?


    '''
    Top BC's
    '''
    # free traction: impose on r3 the condition sigma_zz=0
    A[-1-C[-1]['N'],:] = 0
    B[-1-C[-1]['N'],:] = 0
    A[-1-C[-1]['N'],-1-3*C[-1]['N']:-1-2*C[-1]['N']] = -C[-1]['lambdmu'][-1]*C[-1]['D'][-1,:]
    B[-1-C[-1]['N'],-1-3*C[-1]['N']                 ]   = +C[-1]['lambd']  [-1]

    # sigma_xz (R4=0) sigma_xz==0
    #  disp(C[-1].solid)
    #  if (C[-1].solid==1)
    A[-1,:] = 0
    A[-1,-1] = 1
    B[-1,:] = 0
    
    return A,B

def EigRW(C,A,B,Nmode):
    nerr = 0
    Nl = C[0]['Nl']

    ## solve GEP
    # [Ll,Lr] = RaySel(C);

    A[np.isnan(A)] = 0
    A[np.isinf(A)] = 0
    B[np.isnan(B)] = 0
    B[np.isinf(B)] = 0
    Anew=A;
    Bnew=B;
    
    '''
    % Anew = (Ll*A*Lr);

    % Bnew = (Ll*B*Lr);
    % if min(C(1).mu1D)==0 % if there is a water layer at the top
    % %     % impose continuity condition for r2
    % disp('got water')
    %     n=C(end).N;
    %    Anew(end-2*n+1,:)=0;
    %    Anew(end-2*n+1,end-2*n-C(end-1).N)=-1;
    %    Anew(end-2*n+1,end-2*n)=1;
    %    Bnew(end-2*n+1,:)=0;
    %    
    %    % R2=0
    %    Anew(end-C(end).N,:)=0;
    % Bnew(end-C(end).N,:)=0;
    % Anew(end-C(end).N,end-3*C(end).N+1:end-2*C(end).N) = -C(end).lambdamu(end).*C(end).D(end,:); % )
    % Bnew(end-C(end).N,end-3*C(end).N                 ) = +C(end).lambda  (end);
    % 
    % 
    % end
    '''

    [X,k] = np.linalg.eig(Anew,Bnew)
    %[X,k] = eig(A,B)

    k = diag(k)
    ik = np.where(np.real(k)/np.imag(k)>=1E15 & np.real(k) ~= np.inf & np.real(k) > 0 & C[0]['omega']/np.real[k]>0.1) 
    k = k[ik]
    X = X(:,ik)
    [k,ii] = sort(k,'descend')
    X = X(:,ii);

    # figure(301)
    # plot(C(1).omega/2/pi./k,'o');grid on
    # pause

    disp(['# of Rayleigh modes found: ' int2str(length(k))])
    if isempty(k):
        nerr = 1
        return

    #disp('magic')
    #size(X)
    # X = Lr*X;
    #size(X)
    #figure; plot(X(:,1))
    #disp('magic')
    #figure; imagesc(Lr); 

    n1=0;W=[];bigD=0;
    for i in range(C[0]['Nl']):
        n2 = n1 + C(i).N;
        bigD(n1+1:n2,n1+1:n2)=C(i).D;
        [~,W(n1+1:n2)]=clencurt(C(i).N-1);
        W(n1+1:n2) = C(i).H/2*W(n1+1:n2);
        C(i).nn=[n1 n2];
        n1=n2;

    C(1).BigD = bigD;
    C(1).W = W;


    %% eigenfunctions:
    C(1).ux = zeros(C(1).Ntot,1);
    C(1).uz = zeros(C(1).Ntot,1);
    C(1).r3 = zeros(C(1).Ntot,1);
    C(1).r4 = zeros(C(1).Ntot,1);
    C(1).szz= zeros(C(1).Ntot,1);

    for ik=1:min([length(k) Nmode])
        n=0;
        toto1=X(0*C(1).N+1:1*C(1).N,ik)';
        toto2=X(1*C(1).N+1:2*C(1).N,ik)';
        toto3=X(2*C(1).N+1:3*C(1).N,ik)';
        toto4=X(3*C(1).N+1:4*C(1).N,ik)';
        for i=2:Nl
          n = n + 4*C(i-1).N;
          toto1 = [toto1 , X(n+0*C(i).N+1:n+1*C(i).N,ik)'];
          toto2 = [toto2 , X(n+1*C(i).N+1:n+2*C(i).N,ik)'];
          toto3 = [toto3 , X(n+2*C(i).N+1:n+3*C(i).N,ik)'];
          toto4 = [toto4 , X(n+3*C(i).N+1:n+4*C(i).N,ik)'];
        end
        
        normind=C(1).Ntot;
    %     for ii=Nl:-1:1; if (C(1).solid(ii)~=1); normind=normind-C(ii).N; break ; end; end;
        %[null,normind] = min(C.solid1D); normind = normind-1;
        %if (min(C(1).solid1D)==1); normind=length(C(1).solid1D); end;
        %if (max(C(1).solid1D)==0); normind=length(C(1).solid1D); end;
    %     normalize=max(abs(toto2))*sign(toto2(end));%(normind);
    %     if ik==1;normalize=toto2(end);end
    %     if (normalize==0); normalize = max(toto2); end;
        normalize=toto2(end);
        C(1).ux(:,ik) = toto1/normalize;
        C(1).uz(:,ik) = toto2/normalize;
        C(1).r3(:,ik) = toto3/normalize;
        C(1).r4(:,ik) = toto4/normalize;
    '''
    %     if ik==1&& min(toto2)*max(toto2)<0&&min(toto2)<-0.1;
    %             disp('modes flipped')
    %     end
    %     CC.ux(:,ik) = toto1;
    %     CC.uz(:,ik) = toto2;
    %     CC.r3(:,ik) = toto3./CC.ux(end,ik);
    %     CC.r4(:,ik) = toto4./CC.ux(end,ik);
    %     CC.uz(:,ik) = CC.uz(:,ik)./CC.ux(end,ik);CC.ux(:,ik)=CC.ux(:,ik)./CC.ux(end,ik); 
    '''
        C(1).kr(ik)=k(ik);

        C(1).szz(:,ik) = 4*(C(1).mu1D.*(C(1).lambda1D+C(1).mu1D)./C(1).lambdamu1D).*(C(1).BigD*C(1).uz(:,ik))' ...
            + C(1).lambda1D./C(1).lambdamu1D.*C(1).r3(:,ik)';
        
        
    #end
        return C, nerr
