/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "sineWaveToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sineWaveToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, sineWaveToCell, word);
    addToRunTimeSelectionTable(topoSetSource, sineWaveToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::sineWaveToCell::usage_
(
    sineWaveToCell::typeName,
    "\n    Usage: sineWaveToCell (p1X p1Y p1Z) (p2X p2Y p2Z) radius\n\n"
    "    Select all cells with cell centre within bounding triangle\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sineWaveToCell::combine(topoSet& set, const bool add) const
{
    // const vector axis = p2_ - p1_;
    // const scalar magAxis2 = magSqr(axis);

    const pointField& ctrs = mesh_.cellCentres();

    forAll(ctrs, celli)
    {
        
        scalar yloc = ctrs[celli].component(1);
        scalar waveHeight = offset_+amplitude_*sin(((2*yloc*3.14159265359)/lambda_)+phase_);
        scalar xloc = ctrs[celli].component(0);
        if (xloc < waveHeight)
        {
            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sineWaveToCell::sineWaveToCell
(
    const polyMesh& mesh,
    const scalar offset,
    const scalar amplitude,
    const scalar lambda,
    const scalar phase
)
:
    topoSetSource(mesh),
    offset_(offset),
    amplitude_(amplitude),
    lambda_(lambda),
    phase_(phase)
    
{}


Foam::sineWaveToCell::sineWaveToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    offset_(readScalar(dict.lookup("offset"))),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    lambda_(readScalar(dict.lookup("lambda"))),
    phase_(readScalar(dict.lookup("phase")))
{}


Foam::sineWaveToCell::sineWaveToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    offset_(readScalar(checkIs(is))),
    amplitude_(readScalar(checkIs(is))),
    lambda_(readScalar(checkIs(is))),
    phase_(readScalar(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sineWaveToCell::~sineWaveToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sineWaveToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
     {
    //     Info<< "    Adding cells with centre within triangle, with p1 = "
    //         << p1_ << ", p2 = " << p2_ << " and radius1 = " << radius1_ 
    //         << ", radius2 = " << radius2_ << endl;

         combine(set, true);
    }
     else if (action == topoSetSource::DELETE)
     {
    //     Info<< "    Removing cells with centre within triangle, with p1 = "
    //         << p1_ << ", p2 = " << p2_ << " and radius1 = " << radius1_ 
    //         << ", radius2 = " << radius2_ << endl;

         combine(set, false);
     }
}


// ************************************************************************* //
