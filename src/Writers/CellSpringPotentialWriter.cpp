/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "CellSpringPotentialWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellSpringPotentialWriter<ELEMENT_DIM, SPACE_DIM>::CellSpringPotentialWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellSpringPotential.dat")
{
    this->mVtkCellDataName = "Cell Spring Potential";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellSpringPotentialWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    // double cell_id = pCell->GetCellId();
    double cell_spring_potential = pCell->GetCellData()->GetItem("cell_spring_potential");
    return cell_spring_potential;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellSpringPotentialWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{       
    unsigned cell_id = pCell->GetCellId();
    double cell_spring_potential = pCell->GetCellData()->GetItem("cell_spring_potential");

    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    *this->mpOutStream << " " << cell_id << " " << location_index << " " << std::fixed << std::setprecision(12) << cell_spring_potential;


}

// Explicit instantiation
template class CellSpringPotentialWriter<1,1>;
template class CellSpringPotentialWriter<1,2>;
template class CellSpringPotentialWriter<2,2>;
template class CellSpringPotentialWriter<1,3>;
template class CellSpringPotentialWriter<2,3>;
template class CellSpringPotentialWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellSpringPotentialWriter)
