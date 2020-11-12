// Haematocrit transport equation


gammaDot = pow(2,0.5)*mag(symm(fvc::grad(U)));
// that is the usual definition of gammaDot

sourceC         = fvc::laplacian(linKc*Kc*sqr(a)*sqr(H), gammaDot);

sourceCnonlin   = fvc::laplacian(nonlinKc*nonlinKcRatio*Kmu*H*sqr(a)*sqr(H),
                                 gammaDot);

sourceV         = fvc::laplacian(Kmu*sqr(a)
                                 * gammaDot*sqr(H)/laminarTransport.nu(),
                                 laminarTransport.nu());

// The following are calculated using the explicit method here 
// and written out as fields for debugging, but the actual H equation 
// uses fvm:: for implicit solution:

diffH           = fvc::laplacian(linKc*Kc*sqr(a)*H*gammaDot, H);

diffHnonlin     = fvc::laplacian(nonlinKc*nonlinKcRatio*Kmu*H*sqr(a)*H*gammaDot, H);

fvScalarMatrix HEqn
(
 fvm::ddt(H)
 + fvm::div(phi, H)
 - fvm::laplacian(linKc*Kc*sqr(a)*H*gammaDot, H) // diffH
 - fvm::laplacian(nonlinKc*nonlinKcRatio*Kmu*H*sqr(a)*H*gammaDot, H) // diffHnonlin
 ==
 sourceC 
 + sourceCnonlin 
 + sourceV
 );
