#ifndef TEST3DBOXMODEL_HPP_
#define TEST3DBOXMODEL_HPP_

#include <cxxtest/TestSuite.h>
#include "Timer.hpp"

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "TrianglesMeshReader.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
// #include "SimpleWntCellCycleModel.hpp"
#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StromalCellMutationState.hpp"
#include "CellsGenerator.hpp"
#include "RandomMotionForce.hpp"
#include "PeriodicCryptModelInteractionForceWithGhostNodes.hpp"
// #include "PeriodicBendingForce3dHeightWithGhostNodes.hpp"
// #include "SloughingCellKiller3DWithGhostNodes.hpp"
// #include "AnoikisCellKiller3DWithGhostNodes.hpp"
// #include "UniformCellKiller3dWithGhostNodes.hpp"
#include "PeriodicBoxBoundaryCondition3dGhosts.hpp"
#include "PeriodicStromalBoxBoundaryCondition3d.hpp"
#include "TrianglesMeshWriter.hpp"
#include "Debug.hpp"

#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellAncestorWriter.hpp"
// #include "CellPopulationEpithelialWriter.hpp"
#include "VoronoiDataWriter.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "ForwardEulerNumericalMethod.hpp"
// #include "PlanarDivisionRule.hpp"
#include "DriftPreventForce.hpp"

#include "PeriodicNeighbourModifier.hpp"
// #include "PeriodicNeighbourNonGhostModifier.hpp"

// #include "MeshModifier.hpp"
// #include "MeshRemeshModifier.hpp"
// #include "PeriodicRemeshCellsModifier.hpp"

#include "CellSpringPotentialWriter.hpp"

// Tests have been copied directly from TestOffLatticeSimulation3d.hpp - most likely that bits will get changed
// and broken by me along the way

class Test3dBoxModel : public AbstractCellBasedTestSuite
{
private:

public:
    
    /* A cube that consists of a block of stromal cells, and a single layer of epithelial cells.
     * Target curvature for the basement membrane is zero, and multiple division events occur.
     */
    void TestPeriodicCubeWithGhosts() throw (Exception)
    {
    	RandomNumberGenerator::Instance()->Reseed(1);

        std::vector<Node<3>*> nodes;

        std::string output_directory = "Packing/Random_Static_SCC/Noise_0.01";

        unsigned width = 10;	   // x
        unsigned height = 10;      // y
        unsigned ghosts_bottom = 1;       // ghosts > depth
        unsigned ghosts_top = 2;       // ghosts > depth
        unsigned num_tissue_depth = 5;
        unsigned depth = num_tissue_depth + (ghosts_bottom + ghosts_top) + 1;        // z

        double rand_noise = 0.01;
        double spring_max_length = 1.5;
        double rand_spread = 1;

        // Initialise the tissue in an equilibrum state
        double width_space = 1.0;
        double height_space = 1.0;
        double ghost_sep = 1.0;
        double depth_space = 1.0; //Magic number for z-spaceing... 
        // double depth_space = 1.0; //Magic number for z-spaceing... 
        unsigned cell_iter = 0;

        double periodic_width = (double) (width+0.0)*width_space;
        double periodic_height = (double) (height+0.0)*height_space;

        double tissue_base = 5.0; //Hieght of tissue to prevent drift

        bool isGhost = false;
        double x_coordinate, y_coordinate, z_coordinate;

        double spring_strength = 20.0;

        double time_step = 0.001;
        double end_time = 4;
        double plot_step = 10.0;
        
        std::vector<unsigned> ghost_node_indices, real_node_indices;
        for (unsigned k=0; k<(depth + 1 ); k++) //+1 puts a layer of ghosts on the bottom
        {
            isGhost = false;
            if(k < ghosts_bottom || k > ghosts_bottom + num_tissue_depth)
            {
                isGhost = true;
            }
            if(k == depth)
            {
                isGhost = true;
            }

            for (unsigned j=0; j<height; j++)
            {    
                for (unsigned i=0; i<width; i++)
                {
                        
                    c_vector<double, 3> node_i_new_location;
                    
                    x_coordinate = (double) ( i + rand_spread*(2.0*RandomNumberGenerator::Instance()->ranf()-1.0 ))*width_space;
                    y_coordinate = (double) ( j + rand_spread*(2.0*RandomNumberGenerator::Instance()->ranf()-1.0 ))*height_space;
                
                    // x_coordinate = (double) ( (i + 0.5*(j%2 + k%2))  )*width_space;
                    // y_coordinate = (double) ( j                      )*height_space;
                
                  

                    if( k == depth)
                    {
                        z_coordinate = (double) tissue_base + ((-1.0) + rand_spread*(2.0*RandomNumberGenerator::Instance()->ranf()-1.0))*depth_space;
                        // z_coordinate = (double) tissue_base + (-1.0)*depth_space;

                    }
                    else
                    {
                        z_coordinate = (double) tissue_base + (k      + rand_spread*(2.0*RandomNumberGenerator::Instance()->ranf()-1.0))*depth_space;
                        // z_coordinate = (double) tissue_base + (k )*depth_space;
                    }    
                    if( pow(x_coordinate - 0.5*periodic_width,2)+ pow(y_coordinate - 0.5*periodic_height ,2) <= pow(1.0,2) )
                    {
                    }
                        
                    if(isGhost)
                    {
                        // if(i == 0)
                        // {
                        //     x_coordinate = x_coordinate + 0.2;
                        // }
                        // else if(i == width - 1)
                        // {
                        //     x_coordinate = x_coordinate - 0.2;
                        // }
                        // if(j == 0)
                        // {
                        //     y_coordinate = y_coordinate + 0.2;
                        // }
                        // else if(j == height - 1)
                        // {
                        //     y_coordinate = y_coordinate - 0.2;
                        // }
                        ghost_node_indices.push_back(cell_iter);
                    }
                    else
                    {
                        real_node_indices.push_back(cell_iter);
                    }
                    nodes.push_back(new Node<3>(cell_iter,  false,  x_coordinate, y_coordinate, z_coordinate));
                    cell_iter++;



                    }
                }

        }

        MutableMesh<3,3> mesh(nodes);

        std::vector<CellPtr> cells;

    	// First we sort the real node indices into increasing order (the last ones will correspond to the
    	// epithelial nodes)
    	sort(real_node_indices.begin(), real_node_indices.end());

        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_stromal_state(new StromalCellMutationState);

        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        // Initialise Tissue cells (Stromal)
		for (unsigned i=0; i<real_node_indices.size(); i++)
		{
			//StochasticDurationGenerationBasedCellCycleModel* p_model = new StochasticDurationGenerationBasedCellCycleModel();
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

            CellPtr p_differentiated_cell(new Cell(p_stromal_state, p_model));
            p_differentiated_cell->SetCellProliferativeType(p_differentiated_type);
            p_model->SetDimension(3);

            double birth_time = - 10.0;
            p_differentiated_cell->SetBirthTime(birth_time);

			cells.push_back(p_differentiated_cell);
        }

        DomMeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices); //ghost_sep
        assert(cell_population.GetNumRealCells() != 0);

        OffLatticeSimulation<3> simulator(cell_population);
        
        // Pass an adaptive numerical method to the simulation
        boost::shared_ptr<AbstractNumericalMethod<3,3> > p_method(new ForwardEulerNumericalMethod<3,3>());
        p_method->SetUseAdaptiveTimestep(false);
        simulator.SetNumericalMethod(p_method);


        // Make sure we have a Voronoi tessellation to begin with
        cell_population.CreateVoronoiTessellation();
        
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        // cell_population.AddPopulationWriter<CellPopulationEpithelialWriter>();
        cell_population.AddCellWriter<CellSpringPotentialWriter>();


        //To fix paraview
        cell_population.SetWriteVtkAsPointsDom(true);
        cell_population.SetOutputMeshInVtkDom(true);



        std::map<Node<3>*, c_vector<double,3> > node_locations_before;
        for (std::list<CellPtr>::iterator cell_iter = cell_population.rGetCells().begin();
             cell_iter != cell_population.rGetCells().end();
             ++cell_iter)
        {
            Node<3>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            node_locations_before[p_node] = p_node->rGetLocation();
            //PRINT_VECTOR(node_locations_before[p_node]);
        }

        MAKE_PTR_ARGS(PeriodicBoxBoundaryCondition3dGhosts, boundary_condition, (&cell_population));
        boundary_condition->SetCellPopulationWidth(periodic_width);
        boundary_condition->SetCellPopulationDepth(periodic_height);
        boundary_condition->SetMaxHeightForPinnedCells(0.0);       			      
        boundary_condition->ImposeBoundaryCondition(node_locations_before);
        simulator.AddCellPopulationBoundaryCondition(boundary_condition);

        MAKE_PTR_ARGS(PeriodicStromalBoxBoundaryCondition3d, stromal_boundary_condition, (&cell_population));
        stromal_boundary_condition->SetCellPopulationWidth(periodic_width);
        stromal_boundary_condition->SetCellPopulationDepth(periodic_height);
        stromal_boundary_condition->SetMaxHeightForPinnedCells(0.0);
        stromal_boundary_condition->ImposeBoundaryCondition(node_locations_before);
        simulator.AddCellPopulationBoundaryCondition(stromal_boundary_condition);

		// Create periodic spring force law
        MAKE_PTR(PeriodicCryptModelInteractionForceWithGhostNodes<3>, periodic_spring_force);
        periodic_spring_force->SetUseOneWaySprings(false); //turning this on makes the stromal cells act as ghosts..
        periodic_spring_force->SetCutOffLength(spring_max_length);
        //                     SetEpithelialStromalCellDependentSprings(ind , Ep-Ep, Str-Str, Ep-Str, apcTwoHitStromalMultiplier);
        periodic_spring_force->SetEpithelialStromalCellDependentSprings(true, 1.0,     1.0,     1.0,    1.0);
        periodic_spring_force->SetPeriodicDomainWidth(periodic_width);
        periodic_spring_force->SetPeriodicDomainDepth(periodic_height);
        periodic_spring_force->SetMeinekeSpringStiffness(spring_strength);
        simulator.AddForce(periodic_spring_force);

        // Prevents getting stuck in a local minimums -> used to help break symmetry in cell anoikus
        if(rand_noise > 0.0)
        {    
            MAKE_PTR(RandomMotionForce<3>, p_random_force);
            p_random_force->SetMovementParameter(rand_noise); //0.1 causes dissasociation, 0.001 is not enough
            simulator.AddForce(p_random_force);
        }

        MAKE_PTR(PeriodicNeighbourModifier<3>, p_n_modifier);
        p_n_modifier->SetOutputDirectory(output_directory + "/results_from_time_0");
        p_n_modifier->SetWidth(periodic_width);
        p_n_modifier->SetDepth(periodic_height);
        simulator.AddSimulationModifier(p_n_modifier);

        simulator.SetOutputDirectory(output_directory);	 

        simulator.SetSamplingTimestepMultiple(plot_step);			// Every hour
		simulator.SetEndTime(end_time);
        simulator.SetDt(time_step);

        Timer::Reset();
        simulator.Solve();
        Timer::Print("Time Ellapsed");
        
    }     

};

#endif /*TEST3DBOXMODEL_HPP_*/

