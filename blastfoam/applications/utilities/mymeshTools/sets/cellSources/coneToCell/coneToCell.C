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

#include "coneToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coneToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, coneToCell, word);
    addToRunTimeSelectionTable(topoSetSource, coneToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::coneToCell::usage_
(
    coneToCell::typeName,
    "\n    Usage: coneToCell (p1X p1Y p1Z) (p2X p2Y p2Z) radius\n\n"
    "    Select all cells with cell centre within bounding triangle\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coneToCell::combine(topoSet& set, const bool add) const
{
    const vector axis = p2_ - p1_;
    const scalar magAxis2 = magSqr(axis);

    const pointField& ctrs = mesh_.cellCentres();

    forAll(ctrs, celli)
    {
        vector d = ctrs[celli] - p1_;
        scalar magD = d & axis;

        if ((magD > 0) && (magD < magAxis2))
        {
            scalar rad2 = sqr(magD/magAxis2*radius2_ + (1.-magD/magAxis2)*radius1_);
            scalar d2 = (d & d) - sqr(magD)/magAxis2;
            if (d2 < rad2)
            {
                addOrDelete(set, celli, add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coneToCell::coneToCell
(
    const polyMesh& mesh,
    const vector& p1,
    const vector& p2,
    const scalar radius1,
    const scalar radius2
)
:
    topoSetSource(mesh),
    p1_(p1),
    p2_(p2),
    radius1_(radius1),
    radius2_(radius2)
{}


Foam::coneToCell::coneToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    p1_(dict.lookup("p1")),
    p2_(dict.lookup("p2")),
    radius1_(readScalar(dict.lookup("radius1"))),
    radius2_(readScalar(dict.lookup("radius2")))
{}


Foam::coneToCell::coneToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    p1_(checkIs(is)),
    p2_(checkIs(is)),
    radius1_(readScalar(checkIs(is))),
    radius2_(readScalar(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coneToCell::~coneToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coneToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells with centre within triangle, with p1 = "
            << p1_ << ", p2 = " << p2_ << " and radius1 = " << radius1_ 
            << ", radius2 = " << radius2_ << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells with centre within triangle, with p1 = "
            << p1_ << ", p2 = " << p2_ << " and radius1 = " << radius1_ 
            << ", radius2 = " << radius2_ << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
