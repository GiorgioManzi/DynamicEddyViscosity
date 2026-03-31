/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kOmegaDynamic.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void kOmegaDynamic<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = (k_/(omega_+low_omega_))*(dynamicCmu_/betaStar_);
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> kOmegaDynamic<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>(
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> kOmegaDynamic<BasicMomentumTransportModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>(
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kOmegaDynamic<BasicMomentumTransportModel>::kOmegaDynamic
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            this->coeffDict_,
            0.075
        )
    ),

    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),

    gamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.55
        )
    ),

    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            0.5
        )
    ),

    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.5
        )
    ),

    Cmu_0_ 
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu_0",
            this->coeffDict_,
            0.09
        )
    ),

    dynamicCmuMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "dynamicCmuMin",
            this->coeffDict_,
            0.0
        )
    ),

    dynamicCmuMax_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "dynamicCmuMax",
            this->coeffDict_,
            0.2
        )
    ),

    zeroDecay_
    (
        this->coeffDict_.template lookupOrDefault<Switch>("Zero-Decay", true)
    ),

    physicalProperties_
    (
        IOobject
        (
            "physicalProperties",         
            this->runTime_.constant(),    
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    nStart_
    (
        this->coeffDict().template lookupOrDefault<label>("nStart", 10)
    ),

    window_(
        this->coeffDict().lookupOrDefault("nWindow", 2)
    ),

    low_omega_
    (
        "low_omega",
        dimensionSet(0, 0, -1, 0, 0),
        SMALL
    ),

    low_k_
    (
        "low_k",
        dimensionSet(0, 2, -2, 0, 0),
        SMALL
    ),

    low_Z_
    (
        "low_Z",
        dimensionSet(0, 4, -4, 0, 0),
        SMALL
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    dynamicCmu_
    (
        IOobject
        (
            "dynamicCmu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    CmuStar_
    (
        IOobject
        (
            "CmuStar",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("one", dimless, 1.0)
    ),

    S 
    (
        IOobject
        (
            "S",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor(dimensionSet(0, 0, -1, 0, 0), Zero)
    ),

    gU 
    (
        IOobject
        (
            "gU",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedTensor(dimensionSet(0, 0, -1, 0, 0), Zero)
    ),

    ringIndex_(0)

{
    Info << "kOmegaDynamic constructor called!" << nl;
    Info << "Zero-Decay: " << zeroDecay_ << nl;

////////////////////////////////    PtrLists    ////////////////////////////////////
    
    kPrev_.setSize(window_ - 1);
    omegaPrev_.setSize(window_ - 1);
    UPrev_.setSize(window_ - 1);
    SPrev_.setSize(window_ - 1);
    kprint_.setSize(window_ - 1);
    omegaprint_.setSize(window_ - 1);
    Uprint_.setSize(window_ - 1);
    gUprint_.setSize(window_ - 1);

    for (int i = 0; i < window_ - 1; ++i)
    {
        kPrev_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kPrev_" + name(i),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0,2,-2,0,0), 0)
            )
        );

        omegaPrev_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "omegaPrev_" + name(i),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0,0,-1,0,0), SMALL)
            )
        );

        UPrev_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "UPrev_" + name(i),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedVector(dimensionSet(0,1,-1,0,0), Zero)
            )
        );

        SPrev_.set
        (
            i,
            new volSymmTensorField
            (
                IOobject
                (
                    "SPrev_" + name(i),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedSymmTensor(dimensionSet(0,0,-1,0,0), Zero)
            )
        );

        kprint_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kprint_" + name(i+1),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0,2,-2,0,0), 0)
            )
        );

        omegaprint_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "omegaprint_" + name(i+1),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0,0,-1,0,0), 1e-12)
            )
        );

        Uprint_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "Uprint_" + name(i+1),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedVector(dimensionSet(0,1,-1,0,0), Zero)
            )
        );

        gUprint_.set
        (
            i,
            new volTensorField
            (
                IOobject
                (
                    "gUprint_" + name(i+1),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedTensor(dimensionSet(0,0,-1,0,0), Zero)
            )
        );
    }

    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kOmegaDynamic<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        beta_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        this->coeffDict().readIfPresent("Zero-Decay", zeroDecay_);
        Cmu_0_.readIfPresent(this->coeffDict());
        dynamicCmuMin_.readIfPresent(this->coeffDict());
        dynamicCmuMax_.readIfPresent(this->coeffDict());
        this->coeffDict().readIfPresent("nStart", nStart_);

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicMomentumTransportModel>
void kOmegaDynamic<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local Ref. 
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;

    volScalarField& nut = this->nut_;

    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G(this->GName(), nut*(tgradU() && dev(twoSymm(tgradU()))));

    // Update omega and G at walls
    omega_.boundaryFieldRef().updateCoeffs();   

    // --- Zero-Decay Switch ---
    if (zeroDecay_)
    {
        // If active, use the zero-decay modification 
        CmuStar_ = dynamicCmu_ / Cmu_0_;
    }
    else
    {
        // If inactive, the factor becomes 1.0 (classical case with free-stream decay)
        CmuStar_ = dimensionedScalar("one", dimless, 1.0);
    }

    // Specific dissipation equation
    tmp<fvScalarMatrix> omegaEqn(
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        gamma_*alpha*rho*G*omega_/(k_ + low_k_)
      - fvm::SuSp((2.0/3.0)*gamma_*alpha*rho*divU, omega_)
      - fvm::Sp(CmuStar_*beta_*alpha*rho*omega_, omega_)
      + omegaSource()
      + fvModels.source(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvConstraints.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvConstraints.constrain(omega_);
    bound(omega_, this->omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn(
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(CmuStar_*betaStar_*alpha*rho*omega_, k_)
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);

//////////////////////////////////  Dynamic procedure  ////////////////////////////////////

    const volTensorField& gradU = tgradU();

    S = symm(gradU);
    
    gU = gradU;

    volSymmTensorField sum_uu = symm(U * U);
    volVectorField sum_u = U;

    volSymmTensorField kES_mean = (2/betaStar_) * (k_/(omega_ + low_omega_)) * S;
    volSymmTensorField sum_S = S;

    volScalarField k_mean = k_;
    volScalarField omega_mean = omega_;

    volScalarField squaredS_mean =
        sqr(S.component(symmTensor::XX))
      + sqr(S.component(symmTensor::YY))
      + sqr(S.component(symmTensor::ZZ));


    tgradU.clear();

    for (int i = 0; i < window_ - 1; ++i)
    {
        sum_uu += symm(UPrev_[i] * UPrev_[i]);
        sum_u += UPrev_[i];

        const volSymmTensorField& S_i = SPrev_[i];

        kES_mean += (2/betaStar_)*(kPrev_[i]/(omegaPrev_[i] + low_omega_)) * S_i;
        sum_S += S_i;

        k_mean += kPrev_[i];
        omega_mean += omegaPrev_[i];

        squaredS_mean +=
            sqr(S_i.component(symmTensor::XX))
          + sqr(S_i.component(symmTensor::YY))
          + sqr(S_i.component(symmTensor::ZZ));
    }

    sum_u = sum_u / window_;
    volSymmTensorField uu_average = sum_uu / window_;
    volSymmTensorField u_average_u_average = symm((sum_u) * (sum_u));

    volSymmTensorField L = uu_average - u_average_u_average;

    kES_mean = kES_mean / window_;
    volSymmTensorField S_average = sum_S / window_;

    k_mean = k_mean / window_;
    volScalarField k_mean_T = k_mean + 0.5*(tr(uu_average) - (sum_u & sum_u)); 

    omega_mean = omega_mean / window_;

    squaredS_mean = squaredS_mean / window_;

    volScalarField SSquared_mean =
        sqr(S_average.component(symmTensor::XX))
      + sqr(S_average.component(symmTensor::YY))
      + sqr(S_average.component(symmTensor::ZZ));

    omega_mean = (1 / (k_mean_T + low_k_)) * ((k_mean * omega_mean) + ((this->nu()) / betaStar_) * (squaredS_mean - SSquared_mean));

    volSymmTensorField Z = kES_mean - (2/betaStar_) * (k_mean_T / (omega_mean + low_omega_)) * S_average;

    dynamicCmu_ = (Z && L) / ((Z && Z) + low_Z_);

    // Starting procedure 
    if (this->runTime_.timeIndex() < nStart_)
    {
        dynamicCmu_ = dimensionedScalar("cmuStart", dimless, 0.09);
    }

    // Clipping
    dynamicCmu_ = max(dynamicCmuMin_, min(dynamicCmu_, dynamicCmuMax_));

    // Update stored fields
    Uprint_[ringIndex_] = UPrev_[ringIndex_];
    kprint_[ringIndex_] = kPrev_[ringIndex_];
    omegaprint_[ringIndex_] = omegaPrev_[ringIndex_];
    gUprint_[ringIndex_] = fvc::grad(UPrev_[ringIndex_]);

    // Insert new values
    UPrev_[ringIndex_] = this->U_;
    kPrev_[ringIndex_] = this->k_;
    omegaPrev_[ringIndex_] = this->omega_;
    SPrev_[ringIndex_] = S;

    // Advance index
    ringIndex_ = (ringIndex_ + 1) % (window_ - 1);

    correctNut();

    Info << "Time = " << this->runTime_.timeName() << nl << endl;
}

////////////////////////////////////////////////////////////////////////////////

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
