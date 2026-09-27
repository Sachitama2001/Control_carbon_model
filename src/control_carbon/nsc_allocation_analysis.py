"""Scalar equilibrium reduction and one-sided boundary audit for V4/V5.

For a constant leaf allocation c, C_i=y c_i V/l_i at equilibrium.
At fixed V, NSC balances give a unique N_L through the mobilization identity;
N_S,N_R are positive roots of quadratics. The remaining leaf NSC balance is
one scalar F(V)=0. Roots must satisfy the original branch inequalities.
"""
import numpy as np
from scipy.optimize import brentq

from .nsc_allocation import evaluate_allocation, branch_conditions
from .temperature_nsc_visit import optimum_leaf_carbon, size_dependent_respiration, visit_physiology


def exclusion_bounds(temperature,base,visit,env,reference_days=1.0):
    """Necessary conditions from G=V+paid+NSC turnover at equilibrium.

    Strict failures exclude positive branch AND sliding equilibria, not cycles.
    No dependence on NSC mobilization rates or nonnegative respiration demand.
    """
    target,_,_=optimum_leaf_carbon(temperature,visit,env)
    phys=visit_physiology(temperature,visit,env)
    psat=phys["psat"]
    if psat<=0:
        slope=cap=0.0
    else:
        b=visit.extinction*phys["lue"]*env.ppfd_top/psat
        factor=psat*env.day_length*(3600*12/1e8)
        slope=factor*(visit.sla*2.2/200)*(1-1/np.sqrt(1+b))
        cap=2*factor/visit.extinction*np.log((1+np.sqrt(1+b))/2)
    loss=base.turnover[0]+base.mu0
    far_gain=base.yield_c*.05*slope/loss
    near_min=loss*target/(visit.alloc_ass*(base.yield_c+loss*reference_days))
    return dict(gpp_origin_slope=float(slope),gpp_upper_bound=float(cap),
                far_gain_ratio=float(far_gain),near_min_mobilization=float(near_min),
                far_excluded=bool(far_gain<1),near_and_surfaces_excluded=bool(cap<near_min),
                positive_equilibria_excluded=bool(far_gain<1 and cap<near_min))


def donor_balance(incoming, loss, demand, half):
    """Nonnegative solution of loss*N + demand*N/(half+N) = incoming."""
    if incoming == 0:
        return 0.0
    b = loss*half+demand-incoming
    disc = np.sqrt(b*b+4*loss*incoming*half)
    return (2*incoming*half/(b+disc) if b >= 0 else (disc-b)/(2*loss))


def reconstruct(total, cf, temperature, base, visit, mobilization, variant, rescue_strength=0):
    loss = np.asarray(base.turnover)+base.mu0
    fractions = np.array([cf, (1-cf)*visit.alloc_abg, (1-cf)*(1-visit.alloc_abg)])
    c = base.yield_c*fractions*total/loss
    if rescue_strength:
        if variant!="V7_V4" or not 0<=rescue_strength<=1:
            raise ValueError("rescue strength only for V7_V4")
        critical=.1*200/(visit.sla*2.2)
        hs=rescue_strength*critical*visit.alloc_abg/mobilization.reference_days
        hr=rescue_strength*critical*(1-visit.alloc_abg)/mobilization.reference_days
        c[1]=donor_balance(base.yield_c*fractions[1]*total,loss[1],hs,10*critical)
        c[2]=donor_balance(base.yield_c*fractions[2]*total,loss[2],hr,5*critical)
        c[0]=(base.yield_c*total-loss[1]*c[1]-loss[2]*c[2])/loss[0]
    demand = (size_dependent_respiration(c, temperature, visit) if variant == "V5"
              else np.asarray(base.respiration)*c*base.q10**((temperature-base.t_ref)/10))
    k = np.asarray(mobilization.rates)
    ls, lr, sl, rl = base.transfer
    def ns(nl):
        return np.array([nl,
            donor_balance(ls*nl, loss[1]+k[1]+sl, demand[1], base.substrate_half[1]),
            donor_balance(lr*nl, loss[2]+k[2]+rl, demand[2], base.substrate_half[2])])
    # Solve in a scaled coordinate to retain relative accuracy near the origin.
    u = brentq(lambda u: np.dot(k, ns(u*total/k[0]))/total-1,
               0, 1, xtol=1e-14)
    n = ns(u*total/k[0])
    soil_rates = np.asarray(base.soil_k)*base.soil_q10**((temperature-base.t_ref)/10)
    inputs = base.retention*np.array([np.dot(loss,c), np.dot(loss,n)])
    x = np.zeros(8)
    x[:3], x[4:7] = c, n
    x[7] = inputs[1]/soil_rates[1]
    x[3] = (inputs[0]+base.humification*inputs[1])/soil_rates[0]
    return x


def finite_jacobian(fun, x):
    j = np.empty((len(x), len(x)))
    for i in range(len(x)):
        h = min(1e-5*max(1.0,x[i]), 0.1*x[i])
        if h <= 0:
            raise ValueError("interior Jacobian required")
        delta = np.eye(len(x))[i]*h
        j[:,i] = (fun(x+delta)-fun(x-delta))/(2*h)
    return j


def sign_change_roots(fun, grid):
    """Finite scan; tangencies/disconnected narrow root pairs are not certified."""
    values = [fun(v) for v in grid]
    roots = []
    for a,b,fa,fb in zip(grid[:-1],grid[1:],values[:-1],values[1:]):
        if fa*fb < 0:
            roots.append(brentq(fun,a,b,xtol=1e-14))
        elif fa == 0:
            roots.append(float(a))
    if values[-1] == 0:
        roots.append(float(grid[-1]))
    return roots


def frozen_audit(temperature, base, visit, env, mobilization, variant, points=600):
    if visit.alloc_ass <= .05:
        raise ValueError("three-constant-branch reduction requires alloc_ass>0.05")
    target, _, status = optimum_leaf_carbon(temperature,visit,env)
    if not np.isfinite(target) or target <= 0:
        return dict(temperature=temperature, variant=variant, status=status, roots=[], boundaries=[])
    def field(x,cf):
        return evaluate_allocation(x,temperature,base,visit,env,mobilization,variant,cf)[0]
    def state(v,cf):
        return reconstruct(v,cf,temperature,base,visit,mobilization,variant)
    def residual(v,cf):
        return field(state(v,cf),cf)[4]/v
    # Reference upper bound is generous; reported scan range, not a proof.
    grid = np.geomspace(1e-12, 10.0, points)
    roots = []
    for name,cf in (("far",.05),("near",visit.alloc_ass)):
        for v in sign_change_roots(lambda v:residual(v,cf),grid):
            x=state(v,cf)
            h0,h1=branch_conditions(x,target,mobilization,visit)
            admissible = (h1<0 if name=="far" else h0<0 and h1>0)
            eig=np.linalg.eigvals(finite_jacobian(lambda z:field(z,cf),x))
            roots.append(dict(branch=name, cf=cf, mobilization=v, state=x.tolist(),
                h0=h0,h1=h1,admissible=bool(admissible),
                max_real=float(eig.real.max()), residual=float(np.max(np.abs(field(x,cf))))))
    boundaries=[]
    loss_leaf=base.turnover[0]+base.mu0
    k=np.asarray(mobilization.rates)
    for surface,lo,hi in (("target",1e-7,visit.alloc_ass),
                           ("far_near",.05,visit.alloc_ass)):
        def volume(cf):
            if surface=="target":
                return loss_leaf*target/(base.yield_c*cf)
            return target/(visit.alloc_ass*mobilization.reference_days+base.yield_c*cf/loss_leaf)
        for cf in sign_change_roots(lambda cf:residual(volume(cf),cf),np.linspace(lo,hi,points)):
            v=volume(cf); x=state(v,cf)
            normal=np.zeros(8); normal[0]=1
            if surface=="far_near":
                normal[4:7]=visit.alloc_ass*mobilization.reference_days*k
            minus=visit.alloc_ass if surface=="target" else .05
            plus=0.0 if surface=="target" else visit.alloc_ass
            fm,fp=field(x,minus),field(x,plus)
            nm,np_=float(normal@fm),float(normal@fp)
            kind=("attracting_sliding" if nm>0>np_ else
                  "repelling_sliding" if nm<0<np_ else "crossing_or_tangent")
            # Eliminate C_L via h=0; tangent dynamics use seven remaining stocks.
            def tangent_field(z):
                xx=np.r_[0.0,z]
                xx[0]=(target if surface=="target" else
                       target-visit.alloc_ass*mobilization.reference_days*np.dot(k,xx[4:7]))
                f0,f1=field(xx,0.0),field(xx,1.0)
                effective=-np.dot(normal,f0)/np.dot(normal,f1-f0)
                return field(xx,effective)[1:]
            eig=np.linalg.eigvals(finite_jacobian(tangent_field,x[1:]))
            boundaries.append(dict(surface=surface,cf=cf,mobilization=v,state=x.tolist(),
                normal_minus=nm,normal_plus=np_,classification=kind,
                tangent_max_real=float(eig.real.max()),
                residual=float(np.max(np.abs(field(x,cf))))))
    return dict(temperature=temperature,variant=variant,status="valid",roots=roots,
                boundaries=boundaries,scan_total_range=[float(grid[0]),float(grid[-1])],
                scan_points=points)


def bare_invasion(temperature,base,visit,env,mobilization):
    """Far-branch linearization in order C_L,C_S,C_R,N_L,N_S,N_R."""
    cf=.05
    c=np.array([cf,(1-cf)*visit.alloc_abg,(1-cf)*(1-visit.alloc_abg)])
    k=np.asarray(mobilization.rates)
    loss=np.asarray(base.turnover)+base.mu0
    a=np.zeros((6,6))
    a[:3,:3]=-np.diag(loss)
    a[:3,3:]=base.yield_c*np.outer(c,k)
    a[3:,3:]=-np.diag(loss+k)
    for i,j,rate in ((0,1,base.transfer[0]),(0,2,base.transfer[1]),
                     (1,0,base.transfer[2]),(2,0,base.transfer[3])):
        a[i+3,i+3]-=rate; a[j+3,i+3]+=rate
    a[3,0]=exclusion_bounds(temperature,base,visit,env)["gpp_origin_slope"]
    eig=np.linalg.eigvals(a)
    return dict(max_real=float(eig.real.max()),matrix=a.tolist())


def integrate_far_branch(initial,temperature,base,visit,env,mobilization,variant,
                         years=2000,max_step=30,rtol=1e-8,atol=1e-12):
    """Positive log integration stops at far/near surface; never crosses silently."""
    from scipy.integrate import solve_ivp
    if variant not in ("V4","V5"):
        raise ValueError("rescue threshold requires its own event handling")
    x0=np.asarray(initial,float)
    if x0.shape!=(8,) or not np.isfinite(x0).all() or min(x0)<=0:
        raise ValueError("positive finite eight-state initial condition required")
    target,_,_=optimum_leaf_carbon(temperature,visit,env)
    if branch_conditions(x0,target,mobilization,visit)[1]>=0:
        raise ValueError("initial state must be strictly in far branch")
    def fun(t,z):
        x=np.exp(z[:8])
        f,d=evaluate_allocation(x,temperature,base,visit,env,mobilization,variant,.05)
        return np.r_[f/x,d["external"]]
    def boundary(t,z):
        return branch_conditions(np.exp(z[:8]),target,mobilization,visit)[1]
    boundary.terminal=True; boundary.direction=1
    def bare(t,z):
        return np.logaddexp.reduce(z[[0,1,2,4,5,6]])-np.log(1e-12)
    bare.terminal=True; bare.direction=-1
    def domain(t,z):
        return min(z[:8])+600
    domain.terminal=True; domain.direction=-1
    def jac(t,z):
        j=np.zeros((9,9)); h=1e-5
        for i in range(8):
            d=np.eye(9)[i]*h
            j[:,i]=(fun(t,z+d)-fun(t,z-d))/(2*h)
        return j
    sol=solve_ivp(fun,(0,years*365),np.r_[np.log(x0),0],method="Radau",
                  jac=jac,rtol=rtol,atol=atol,max_step=max_step,
                  events=(boundary,bare,domain),dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    times=np.unique(np.r_[np.arange(0,sol.t[-1],365),sol.t[-1]])
    values=sol.sol(times)
    audit_times=np.unique(np.r_[sol.t,times,
        *(sol.t[:-1]+q*np.diff(sol.t) for q in (.25,.5,.75))])
    audit_z=sol.sol(audit_times); audit_x=np.exp(audit_z[:8].T)
    if not np.isfinite(audit_x).all() or np.any(audit_x<=0):
        raise RuntimeError("positivity audit failed")
    err=np.max(np.abs(audit_x.sum(axis=1)-x0.sum()-audit_z[8]))
    return dict(times=times,states=np.exp(values[:8].T),external_integral=values[8],
                minimum_stock=float(audit_x.min()),integrated_budget_error=float(err),
                boundary_event=bool(len(sol.t_events[0])),near_bare_event=bool(len(sol.t_events[1])),
                domain_event=bool(len(sol.t_events[2])))


def rescue_audit(temperature,base,visit,env,mobilization,points=400):
    """V7 on V4: source rescue below critical leaf, far allocation branch only.

    A Filippov boundary candidate is not a native daily fixed point.
    """
    critical=.1*200/(visit.sla*2.2)
    target,_,_=optimum_leaf_carbon(temperature,visit,env)
    def field(x,s):
        return evaluate_allocation(x,temperature,base,visit,env,mobilization,
                                   "V7_V4",.05,s)[0]
    def state(v,s):
        return reconstruct(v,.05,temperature,base,visit,mobilization,"V7_V4",s)
    roots=[]
    for v in sign_change_roots(lambda v:field(state(v,1),1)[4]/v,np.geomspace(1e-12,10,points)):
        x=state(v,1)
        eig=np.linalg.eigvals(finite_jacobian(lambda z:field(z,1),x))
        roots.append(dict(state=x.tolist(),admissible=bool(x[0]<critical and
            branch_conditions(x,target,mobilization,visit)[1]<0),
            max_real=float(eig.real.max()),residual=float(np.max(np.abs(field(x,1))))))
    upper=(base.turnover[0]+base.mu0)*critical/(base.yield_c*.05)
    def constrained(s):
        v=brentq(lambda v:state(v,s)[0]-critical,1e-16,upper*(1+1e-12),xtol=1e-15)
        return state(v,s)
    boundaries=[]
    for s in sign_change_roots(lambda s:field(constrained(s),s)[4],np.linspace(0,1,points)):
        x=constrained(s)
        fm,fp=field(x,1),field(x,0)
        def tangent(z):
            xx=np.r_[critical,z]; f0,f1=field(xx,0),field(xx,1)
            strength=-f0[0]/(f1[0]-f0[0])
            return field(xx,strength)[1:]
        eig=np.linalg.eigvals(finite_jacobian(tangent,x[1:]))
        boundaries.append(dict(state=x.tolist(),effective_rescue=s,
            far_branch_valid=bool(branch_conditions(x,target,mobilization,visit)[1]<0),
            normal_below=float(fm[0]),normal_above=float(fp[0]),
            attracting=bool(fm[0]>0>fp[0]),tangent_max_real=float(eig.real.max()),
            residual=float(np.max(np.abs(field(x,s))))))
    return dict(temperature=temperature,roots=roots,boundaries=boundaries,
                critical_leaf=critical,scope="V7 on inherited V4, far allocation branch")
