//
// Molekel - Molecular Visualization Program
// Copyright (C) 2006, 2007 Swiss National Supercomputing Centre (CSCS)
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301, USA.
//
// $Author$
// $Date$
// $Revision$
//

/// Element table entry.
struct MolekelElement
{
    int atomicNumber;
    const char* symbol;
    double covalentRadius;
    double bondOrderRadius;
    double vdwRadius;
    int maxBondValence;
    double mass;
    double electroNegativity;
    double ionizationPotential;
    double electronicAffinity;
    double red, green, blue;
    const char* name;
};

/// Returns a reference to the element table.
const MolekelElement* GetElementTable();

/// Returns element table size.
int GetElementTableSize();
