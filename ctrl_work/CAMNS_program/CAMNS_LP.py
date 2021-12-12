# Generated with SMOP  0.41-beta
from libsmop import *
# CAMNS_LP.m

    #=====================================================================
# Programmers: 
# Wing-Kin Ma, E-mail: wkma@ieee.org
# Tsung-Han Chan, E-mail: chantsunghan@gmail.com
# Date: Sept. 03, 2009
# -------------------------------------------------------
# Reference: 
# T.-H. Chan, W.-K. Ma, C.-Y. Chi, and Y. Wang, ``A convex analysis 
# framework for blind separation of non-negative sources," 
# IEEE Trans. Signal Process., vol. 56, no. 10, pp. 5120-5134, Oct. 2008.
    
    #======================================================================
# A practical implementation of the CAMNS-LP method
    
    # function hS = CAMNS_LP(X,N)
#======================================================================
# Output: 
# hS is the L-by-N extracted soruce matrix, where L is the data length.
#---------------------------------------------------
# Inputs:
# X is the L-by-M observation matrix, where M is the number of
# observations.
# N is the number of sources. 
#========================================================================
    
    
@function
def CAMNS_LP(X=None,N=None,*args,**kwargs):
    varargin = CAMNS_LP.varargin
    nargin = CAMNS_LP.nargin

    #----------- Define default parameters------------------
    TOL_LP=0.001
# CAMNS_LP.m:29
    
    TOL_EXT=1e-06
# CAMNS_LP.m:30
    
    TOL_ZEROS=1e-06
# CAMNS_LP.m:31
    
    #-----------(Modification 1) Remove All Zero vectors in Data Set X--------
# It is not meaningful to analyze the zero mixtures.
    L,M=size(X,nargout=2)
# CAMNS_LP.m:35
    
    index=find(sum(abs(X.T)) >= TOL_ZEROS)
# CAMNS_LP.m:36
    Xn=X(index,arange())
# CAMNS_LP.m:37
    LL=length(index)
# CAMNS_LP.m:38
    #-----------Affine Set Fitting [Proposition 1]--------------
    d=dot(Xn,ones(M,1)) / M
# CAMNS_LP.m:41
    # For computational efficiency, we use SVD of R instead of EVD of R*R'.
    C,Sigma,V=svds(Xn - dot(d,ones(1,M)),N - 1,'L',nargout=3)
# CAMNS_LP.m:43
    #--------LP Extreme-Point Finding Algorithm [Table 1]---------------
#------------Step 1
    el=0
# CAMNS_LP.m:47
    
    Q1=zeros(LL,1)
# CAMNS_LP.m:48
    hS=[]
# CAMNS_LP.m:49
    lp_cnt=0
# CAMNS_LP.m:50
    
    while el < N:

        #--------------Step 2
        w=randn(LL,1)
# CAMNS_LP.m:53
        w2=dot(Q1.T,w)
# CAMNS_LP.m:53
        r=w - dot(Q1,w2)
# CAMNS_LP.m:54
        b=dot(- C.T,r)
# CAMNS_LP.m:56
        pars.fid = copy(0)
# CAMNS_LP.m:57
        A=- C.T
# CAMNS_LP.m:58
        c=copy(d)
# CAMNS_LP.m:58
        K.l = copy(LL)
# CAMNS_LP.m:58
        tic
        x1,alpha1=sedumi(A,b,c,K,pars,nargout=2)
# CAMNS_LP.m:59
        ttime=copy(toc)
# CAMNS_LP.m:59
        lp_cnt=lp_cnt + 1
# CAMNS_LP.m:59
        fprintf('%dth LP: running time= %2.5f\n',lp_cnt,ttime)
        tic
        x2,alpha2=sedumi(A,- b,c,K,pars,nargout=2)
# CAMNS_LP.m:62
        ttime=copy(toc)
# CAMNS_LP.m:62
        lp_cnt=lp_cnt + 1
# CAMNS_LP.m:62
        fprintf('%dth LP: running time= %2.5f\n',lp_cnt,ttime)
        if el == 0:
            # (Modification 2) To play safe, employ the extreme point
        # validation (Lemma 6) to check the obtained optimal solutions.
        # We will reject the solution if it's not an extreme point (not commonplace by our experience).
            if is_ext_pt(C,d,alpha1,TOL_EXT):
                hS=concat([hS,dot(C,alpha1) + d])
# CAMNS_LP.m:70
                fprintf('Find a new extreme pt (minimizing LP).\n')
            if is_ext_pt(C,d,alpha2,TOL_EXT):
                hS=concat([hS,dot(C,alpha2) + d])
# CAMNS_LP.m:73
                fprintf('Find a new extreme pt (maximizing LP).\n')
        else:
            p_star=abs(dot(r.T,(dot(C,alpha1) + d)))
# CAMNS_LP.m:76
            q_star=abs(dot(r.T,(dot(C,alpha2) + d)))
# CAMNS_LP.m:77
            # the LP solution. Hence p_star or q_star may not be exactly zero,
        # even though they are supposed to be zero theoretically. A threshold is
        # in place to decide the acceptance/rejection of the obtained solutions.
        # Also, extreme point validation is employed just to play safe.
            if p_star / (dot(norm(r),norm(dot(C,alpha1) + d))) >= TOL_LP:
                if is_ext_pt(C,d,alpha1,TOL_EXT):
                    hS=concat([hS,dot(C,alpha1) + d])
# CAMNS_LP.m:85
                    fprintf('Find a new extreme pt (minimizing LP).\n')
            if q_star / (dot(norm(r),norm(dot(C,alpha2) + d))) >= TOL_LP:
                if is_ext_pt(C,d,alpha1,TOL_EXT):
                    hS=concat([hS,dot(C,alpha2) + d])
# CAMNS_LP.m:90
                    fprintf('Find a new extreme pt (maximizing LP).\n')
        el=size(hS,2)
# CAMNS_LP.m:95
        #------------Step 6
        if el > 0:
            Q1,R=qr(hS,0,nargout=2)
# CAMNS_LP.m:97

    
    #--------------when el>N (Modification 4)----------------
# When the number of obtained extreme points happens to be greater than 
# the number of true sources (due to numerical errors or violation of
# assumptions), we truncate by selecting the optimal solution at the 
# last run that yields a higher optimal value
    if el > N:
        hS=dot((p_star > q_star),hS(arange(),arange(1,N,1))) + dot((p_star < q_star),hS(arange(),concat([arange(1,N - 1),N + 1])))
# CAMNS_LP.m:105
    
    Y=zeros(L,N)
# CAMNS_LP.m:107
    Y[index,arange()]=hS
# CAMNS_LP.m:108
    hS=copy(Y)
# CAMNS_LP.m:109