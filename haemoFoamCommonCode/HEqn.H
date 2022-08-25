// Haematocrit transport equation

// Strain Rate Tensor E:
volSymmTensorField E = symm(fvc::grad(U));

// Find the Eigenvectors and Eigenvalues of E:
Ei = eigenVectors(E);
Ev = eigenValues(E);

// We need to find the Eigenvector with the smallest magnitude,
// So we create a helper volVectorField with all positive components
// from the Eigenvalue field
// and then find the minimum:
forAll(Ev,index)
{
Evmag[index] = cmptMag(Ev[index]);
minEv[index] = 0;
for (int n=1; n<3; n=n+1)
{
    if (Evmag[index].component(minEv[index]) > Evmag[index].component(n))
    {
        minEv[index] = n;
    }
};
d3[index].component(0) = Ei[index].component(minEv[index]*3);
d3[index].component(1) = Ei[index].component(minEv[index]*3+1);
d3[index].component(2) = Ei[index].component(minEv[index]*3+2);

if (mag(U[index]) == 0)
{
    d1[index] = vector(0, 0, 0);
} else {
d1[index] = U[index]/mag(U[index]);
};

d2[index] = d3[index] ^ d1[index];
};

// Evmag.component(0) = mag(Ev.component(0));
// Evmag.component(1) = mag(Ev.component(1));
// Evmag.component(2) = mag(Ev.component(2));

gammaDot = pow(2,0.5)*mag(E);
// that is the usual definition of gammaDot
//
//
tensor Z = tensor(1, 0, 0,   0, 1, 0,   0, 0, 1);

gradGammaDot = fvc::grad(gammaDot);
divGammaDotZ = fvc::div(gammaDot*Z.T());

sourceCorig          = fvc::laplacian(linKc*Kc*sqr(a)*sqr(H), gammaDot);
sourceCHgamma        = fvc::laplacian(linKc*Kc*sqr(a)*H, H*gammaDot);
sourceCHgammaGrad    = fvc::div(linKc*Kc*sqr(a)*H*fvc::grad(H*gammaDot));
sourceCHgammaZ       = fvc::div(linKc*Kc*sqr(a)*H*fvc::div(H*gammaDot*Z));
sourceCgammaZ        = fvc::div(linKc*Kc*sqr(a)*sqr(H)*fvc::div(gammaDot*Z));
sourceC              = sourceCorig;

double   diffHactive  = 1.0;

// sourceCnonlin   = fvc::laplacian(nonlinKc*nonlinKcRatio*Kmu*H*sqr(a)*sqr(H),
//                                  gammaDot);

sourceV         = fvc::laplacian(Kmu*sqr(a)
                                 * gammaDot*sqr(H)
                                 * 1/laminarTransport.nu(),
                                 laminarTransport.nu());

// sourceV         = fvc::div(Kmu*sqr(a) * H
//                                  * (gammaDot*H*Z)
//                                  & 1/laminarTransport.nu()
//                                  * fvc::grad(laminarTransport.nu()));

// The following are calculated using the explicit method here
// and written out as fields for debugging, but the actual H equation 
// uses fvm:: for implicit solution:

diffH           = fvc::laplacian(linKc*Kc*sqr(a)*H*gammaDot, H);

// diffHnonlin     = fvc::laplacian(nonlinKc*nonlinKcRatio*Kmu*H*sqr(a)*H*gammaDot, H);

fvScalarMatrix HEqn
(
 fvm::ddt(H)
 + fvm::div(phi, H)
 - fvm::laplacian(diffHactive*linKc*Kc*sqr(a)*H*gammaDot, H) // diffH
 // - fvm::laplacian(nonlinKc*nonlinKcRatio*Kmu*H*sqr(a)*H*gammaDot, H)  // diffHnonlin
 ==
 sourceC
 // + sourceCnonlin
 + sourceV
 );
