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

#include "kEpsDyn.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsDyn<BasicMomentumTransportModel>::f2() const
{
	const volScalarField yStar(pow(this->nu()*epsilon_,scalar(0.25))*y_/this->nu()); 
	
    volScalarField eps(epsilon_);
    tmp<volScalarField> Rt = sqr(k_)/(this->nu()*bound(eps, this->epsilonMin_));	
	
    return
        min((scalar(1)-0.3*exp(-sqr(Rt/6.5)))*sqr(scalar(1)-exp(-yStar/3.1)),scalar(1.0));
}


template<class BasicMomentumTransportModel>
void kEpsDyn<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = dynamicCmu_*sqr(k_)/(epsilon_+low_epsilon_);
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kEpsDyn<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kEpsDyn<BasicMomentumTransportModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kEpsDyn<BasicMomentumTransportModel>::kEpsDyn
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

    /*Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),*/
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.5
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.9
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.4
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.4
        )
    ),
	finestra_(this->coeffDict().lookupOrDefault("nPrecedents", 2)), //FILIPPO
	physicalProperties_
    (
        IOobject
        (
            "physicalProperties",         // file name
            this->runTime_.constant(),       // constant directory, e.g., "constant"
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
		//this->physicalProperties_.lookup("nu")
    ),
	
	/*nu_
	(
	    "nu",
		dimViscosity,
		physicalProperties_.lookup("nu")
	),*/
	low_epsilon_
	(
	    "low_epsilon",
		dimensionSet(0, 2, -3, 0, 0),
		SMALL
	),
	low_k_
	(
	    "low_k",
		dimensionSet(0, 2, -2, 0, 0),
		SMALL
	),
	low_M_
	(
	    "low_M",
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
    dynamicCmu_
    (
        IOobject
        (
            "dynamicCmu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
		//dimensionedScalar(dimensionSet(0, 0, 0, 0, 0), 0.09)
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
	/*Sns
    (
        IOobject
        (
            "Sns",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
		dimensionedTensor(dimensionSet(0, 0, -1, 0, 0), Zero)
    ),*/
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
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
		
	y_(wallDist::New(this->mesh_).y()),
	
	f2Field_
    (
        IOobject
        (
            "f2",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        f2()()
    )
{
////////////////////////FILIPPO////////////////////////
	kPrev_.setSize(finestra_-1);
	epsilonPrev_.setSize(finestra_-1);
	UPrev_.setSize(finestra_-1);
	kprint_.setSize(finestra_-1);
	epsilonprint_.setSize(finestra_-1);
	Uprint_.setSize(finestra_-1);
	gUprint_.setSize(finestra_-1);
	for (int i = 0; i < finestra_-1; i++)
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
                dimensionedScalar(dimensionSet(0, 2, -2, 0, 0), 0)
			)
		);
		epsilonPrev_.set
		(
		    i,
			new volScalarField
			(
			    IOobject
                (
				    "epsilonPrev_" + name(i),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0, 2, -3, 0, 0), SMALL)
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
                dimensionedVector(dimensionSet(0, 1, -1, 0, 0), Zero)
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
                dimensionedScalar(dimensionSet(0, 2, -2, 0, 0), 0)
			)
		);
		epsilonprint_.set
		(
		    i,
			new volScalarField
			(
			    IOobject
                (
				    "epsilonprint_" + name(i+1),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0, 2, -3, 0, 0), 1e-12)
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
                dimensionedVector(dimensionSet(0, 1, -1, 0, 0), Zero)
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
                dimensionedTensor(dimensionSet(0, 0, -1, 0, 0), Zero)
			)
		);
    }
////////////////////////FILIPPO////////////////////////
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);


    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kEpsDyn<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        //Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
		//finestra_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void kEpsDyn<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
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
    //tgradU.clear();


    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha*rho*G*epsilon_/(k_+low_k_)
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha*rho*divU, epsilon_)
      - fvm::Sp(C2_*f2()*alpha*rho*epsilon_/(k_+low_k_), epsilon_)
      + epsilonSource()
      + fvModels.source(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    bound(epsilon_, this->epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G - fvm::SuSp(2.0/3.0*alpha*rho*divU, k_)
      - fvm::Sp(alpha*rho*epsilon_/(k_+low_k_), k_)
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);

    //correctNut();
	
////////////////////////FILIPPO////////////////////////	


////////////////////////Equation check////////////////////////

    /*
    tensor A = tensor(1,2,3,4,5,6,7,8,9);
	tensor B = A.T();
	tensor C = tensor(2,7,4,9,1,3,5,2,7);
	scalar ac = A && C;
	
	vector v = vector(1,2,3);
	vector w = vector(4,5,6);
	tensor VW = v*w;
	scalar s3 = v & w;
	
	scalar s1 = 2;
	scalar s2 = 3;
	tensor D = s1/s2*A;
	
	
	forAll(A, i)
	{
		Info << "Component " << i << "of A is equal to: " << A[i] << nl;
	}
	
	forAll(A, i)
	{
		Info << "Component " << i << "of B is equal to: " << B[i] << nl;
	}
	
	Info << "ac is equal to: " << ac << nl;

    if (ac==201)
	{
		Info << "Hand calculation of ac is 201, so && does the inner product as we want!!" << nl;
	}
	
	forAll(VW, i)
	{
		Info << "Component " << i << "of VW is equal to: " << VW[i] << nl;
	}
	
	Info << "s3 is equal to: " << s3 << nl;
	
	forAll(D, i)
	{
		Info << "Component " << i << "of D is equal to: " << D[i] << nl;
	}
	*/
	
	// I had checked whole the operations and the used operators are ok!!
	




////////////////////////Equation check////////////////////////
	const volTensorField& gradU = tgradU();
	S = symm(gradU);       // https://caefn.com/openfoam/tensor-operations look at this why I wrote S in this way
	//Sns = 0.5*(gradU + gradU.T());
	gU = gradU;
	
	
	volSymmTensorField sum_uu = symm(U * U);  
	volVectorField sum_u = U;
	
	volSymmTensorField kES_mean = 2*sqr(k_)/(epsilon_+low_epsilon_)*S;    
	volSymmTensorField sum_S = S;
	
	volScalarField k_mean = k_;
	volScalarField epsilon_mean = epsilon_;
	
	volScalarField squaredS_mean = sqr(tr(S)) - sqr(tr(S));  //completamente a caso, solo per inizializzarlo;
	
	forAll(S, cellI)
    {
          squaredS_mean[cellI] =
          sqr(S[cellI].xx())
        + sqr(S[cellI].yy())
        + sqr(S[cellI].zz());
    }

	tgradU.clear();
	
	
	for (int i = 0; i < finestra_-1; i++)
    {
        sum_uu += symm(UPrev_[i] * UPrev_[i]);    
		sum_u += UPrev_[i];
		
		volTensorField gradUPrev_ = fvc::grad(UPrev_[i]);
		volSymmTensorField SPrev_ = symm(gradUPrev_);
		
		kES_mean += 2*sqr(kPrev_[i])/(epsilonPrev_[i]+low_epsilon_)*SPrev_;
		sum_S += SPrev_;
		
		k_mean += kPrev_[i];
		epsilon_mean += epsilonPrev_[i];
		
		forAll(SPrev_, cellI)
        {
          squaredS_mean[cellI] +=
          sqr(SPrev_[cellI].xx())
		  + sqr(SPrev_[cellI].yy())
		  + sqr(SPrev_[cellI].zz());
        }
    }
	
	sum_u = sum_u/(finestra_);
	volSymmTensorField uu_average = sum_uu/(finestra_);
	volSymmTensorField u_average_u_average = symm((sum_u) * (sum_u));
		
	volSymmTensorField L = uu_average - u_average_u_average;
	
	kES_mean = kES_mean/(finestra_);
	volSymmTensorField S_average = sum_S/(finestra_);

	
	k_mean = k_mean/(finestra_);
	k_mean = k_mean + 0.5*(tr(uu_average) - (sum_u & sum_u));
	
	epsilon_mean = epsilon_mean/(finestra_);
	
	squaredS_mean = squaredS_mean/(finestra_);
	
	volScalarField SSquared_mean = sqr(tr(S_average))-sqr(tr(S_average));
	
	forAll(S_average, cellI)
    {
          SSquared_mean[cellI] +=
          sqr(S_average[cellI].xx())
		  + sqr(S_average[cellI].yy())
		  + sqr(S_average[cellI].zz());
    }
	Info << "Dimensions of squaredS_mean: " << squaredS_mean.dimensions() << nl;
    Info << "Dimensions of SSquared_mean: " << SSquared_mean.dimensions() << nl;
	//Info << "Dimensions of nu_: " << nu_.dimensions() << nl;
	
	epsilon_mean = epsilon_mean + this->nu()*(squaredS_mean - SSquared_mean);
	
	
	
	volSymmTensorField M = kES_mean - 2*sqr(k_mean)/(epsilon_mean+low_epsilon_)*S_average;
	
	Info << "Dimensions of dynamicCmu_: " << dynamicCmu_.dimensions() << nl;
	
	Info << "Dimensions of M: " << M.dimensions() << nl;
    Info << "Dimensions of L: " << L.dimensions() << nl;
	
	dynamicCmu_ = (M && L) / ((M && M)+low_M_);       
	
	
	Info << "max value of dynamicCmu: " << max(dynamicCmu_) << nl;
	Info << "min value of dynamicCmu: " << min(dynamicCmu_) << nl;
	
	if (this->runTime_.timeIndex() < 10)
	{
		dynamicCmu_ = 0.09;
	}
	
	forAll(dynamicCmu_, cellI)
	{
	    if (dynamicCmu_[cellI] > 0.2)
	    {
	    	dynamicCmu_[cellI] = 0.2;
	    }
	    else if (dynamicCmu_[cellI] < 0)
	    {
	    	dynamicCmu_[cellI] = 0;
        }
	}
    //Updating variables
    
	for (int i = 0; i <= finestra_ - 2; i++)
        {
            Uprint_[i] = UPrev_[i];
            kprint_[i] = kPrev_[i];
			epsilonprint_[i] = epsilonPrev_[i];
			gUprint_[i] = fvc::grad(UPrev_[i]);
        }
	
	for (int i = 0; i < finestra_ - 2; i++)
        {
			bound(kPrev_[i], this->kMin_);
			bound(epsilonPrev_[i], this->epsilonMin_);
			
            UPrev_[i] = UPrev_[i+1];
            kPrev_[i] = kPrev_[i+1];
			epsilonPrev_[i] = epsilonPrev_[i+1];
        }
	UPrev_[finestra_ - 2] = this->U_;
    kPrev_[finestra_ - 2] = this->k_;
	epsilonPrev_[finestra_ - 2] = this->epsilon_;
	
	// Ensure the boundary conditions are consistent
    //UPrev_[finestra_ - 2].correctBoundaryConditions();
    //kPrev_[finestra_ - 2].correctBoundaryConditions();
	//epsilonPrev_[finestra_ - 2].correctBoundaryConditions();
	
    correctNut();
	
// time step number (i.e., time index)
Info << "Time Step Index: " << this->runTime_.timeIndex() << endl;

////////////////////////FILIPPO////////////////////////
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
