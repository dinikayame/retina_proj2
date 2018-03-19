#ifndef RETINA_PROJ2_H_
#define RETINA_PROJ2_H_

#include <vector>
#include <math.h>
#include "biodynamo.h"
#include "math_util.h"
#include "matrix.h"
//#include "substance_initializers.h"

namespace bdm {
  /*
  part 2:
  TODO add substance to this model
  use substance as a way to control the way the cells migrate
  */

/* extend into cell class to define the cell types
*/
  BDM_SIM_OBJECT(MyCell, Cell) {
    BDM_SIM_OBJECT_HEADER(MyCellExt, 1, cellType);

  public:
    MyCellExt() {}
    MyCellExt(const std::array<double, 3>& position) : Base(position) {}

    void SetCellType(int t) { cellType[kIdx] = t; }

    int GetCellType() { return cellType[kIdx]; }
    // This function is used by ParaView for coloring the cells by their type
    int* GetCellTypePtr() { return cellType.data(); }

  private:
    vec<int> cellType;
  };

  /* List the extracellular substances

  Part 2: substance will lead cell migration
  use various concentration to lead the cell to meet

  Q: use different substances for each layer?

  substance behavior is stated in a biology module

  */
 enum Substances { kSubstance };

  /*
  new biology module to aid cell migration -> use concentration levels to push
  cells (i.e. diffusion grad)
  Chemotaxis is to move the cells based on the substance stimuli
  */
  struct Chemotaxis : public BaseBiologyModule {
    Chemotaxis() : BaseBiologyModule(gAllBmEvents) {}

    template <typename T>
    void Run(T* cell) {
      static auto dg = GetDiffusionGrid(kSubstance);
      //TODO change value of threshold later
      dg->SetConcentrationThreshold(1e15);

      auto& position = cell->GetPosition();
      std::array<double, 3> gradient;
      dg->GetGradient(position, &gradient);
      //TODO change values if needed
      gradient[0] *= 0.5;
      gradient[1] *= 0.5;
      gradient[2] *= 0.5;

      cell->UpdatePosition(gradient);
    }
  };

  /*
  this is to assign 1 cell to secret the substance artificially at 1 location
  use this to link the cells within the different layers
  i.e. like how the inner plexiform layer has ganglion cells + amacrine cells
  */
  struct kSecretion : public BaseBiologyModule {
    kSecretion() : BaseBiologyModule() {}

    template <typename T>
    void Run(T* cell){
      static auto dg = GetDiffusionGrid(kSubstance);
      //TODO may or may not need this
      //needs to be declared/ used somewhere then uncomment this
      //array<double, 3> secretion_position = {50, 50, 50};
      //to change value
      double amount = 4;
      //increase conc --> cell pos changes by that amount
      dg->IncreaseConcentrationBy(cell->GetPosition(), amount);
    }
  };



    /*
    Ganglion layer cells -> focus on 3 types:
    1. Midget cell -> input from BIPOLAR CELLS + small size
    2. Parasol cell -> input from rod & cones; part of M pathway + large size (dentritic)
    3. Bistratified cell -> input from BIPOLAR & AMACRINE CELLS + medium size
    */
    struct midgetCell : public BaseBiologyModule {
      midgetCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //what will make the cell grow
        if (cell->GetCellType() > 0) {
          array<double, 3> direction {0, 1, 0};
          //make it grow along the y axis
          cell->UpdateMassLocation(direction);
          //to ensure that there is no traction between cells
          cell->SetTractorForce({0,0,0});
          cell->SetCellType(cell->GetCellType()-1);
        }
      }
    };

    struct parasolCell : public BaseBiologyModule {
      parasolCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //what will make the cell grow
        if (cell->GetCellType() > 0) {
          array<double, 3> direction {0, 1, 0};
          //make it grow along the y axis
          cell->UpdateMassLocation(direction);
          //to ensure that there is no traction between cells
          cell->SetTractorForce({0,0,0});
          cell->SetCellType(cell->GetCellType()-1);
        }
      }
    };

    struct bistratifiedCell : public BaseBiologyModule {
      bistratifiedCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //what will make the cell grow
        if (cell->GetCellType() > 0) {
          array<double, 3> direction {0, 1, 0};
          //make it grow along the y axis
          cell->UpdateMassLocation(direction);
          //to ensure that there is no traction between cells
          cell->SetTractorForce({0,0,0});
          cell->SetCellType(cell->GetCellType()-1);
        }
      }
    };

    /*
    Amacrine cells that are in the INNER LIMITING LAYER
    link bipolar and ganglion cells
    cell takes input from ganglion cell to bipolar cell
    So for ganglion cells --> link to the MIDGET & BISTRATIFIED CELLS
    Cells work laterally
    */
    struct amacrineCell : public BaseBiologyModule {
      amacrineCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //what will make the cell grow
        if (cell->GetCellType() > 0) {
          array<double, 3> direction {0, 1, 0};
          //make it grow along the y axis
          cell->UpdateMassLocation(direction);
          //to ensure that there is no traction between cells
          cell->SetTractorForce({0,0,0});
          cell->SetCellType(cell->GetCellType()-1);
        }
      }
    };

    /*bipolar cells
    exists between photoreceptors and ganglion cells
    so it will either synpase with:
    1. photoreceptors -> rods/ cones (i.e. Parasol cells)
    2. horizontal cells
    */
    struct bipolarCell : public BaseBiologyModule {
      bipolarCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //what will make the cell grow
        if (cell->GetCellType() > 0) {
          array<double, 3> direction {0, 1, 0};
          //make it grow along the y axis
          cell->UpdateMassLocation(direction);
          //to ensure that there is no traction between cells
          cell->SetTractorForce({0,0,0});
          cell->SetCellType(cell->GetCellType()-1);
        }
      }
    };

    /*horizontal cells
    part 2:
    cells work laterally
    connected to outputs from rods and cones
    */
    struct horizontalCell : public BaseBiologyModule {
      horizontalCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //what will make the cell grow
        if (cell->GetCellType() > 0) {
          array<double, 3> direction {0, 1, 0};
          //make it grow along the y axis
          cell->UpdateMassLocation(direction);
          //to ensure that there is no traction between cells
          cell->SetTractorForce({0,0,0});
          cell->SetCellType(cell->GetCellType()-1);
        }
      }
    };

    /*cones
    part 2:
    cells are long but rounder and wider than rods
    connect to horizontal cells + bipolar cells + Parasol cells? <check this>
    */
    struct coneCell : public BaseBiologyModule {
      coneCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //what will make the cell grow
        if (cell->GetCellType() > 0) {
          array<double, 3> direction {0, 1, 0};
          //make it grow along the y axis
          cell->UpdateMassLocation(direction);
          //to ensure that there is no traction between cells
          cell->SetTractorForce({0,0,0});
          cell->SetCellType(cell->GetCellType()-1);
        }
      }
    };

    /*rods
    part 2:
    cells are longish -> outer seg + inner seg + nucleus
    connect to horizontal cells + bipolar cells + Parasol cells? <check this>
    */
    struct rodCell : public BaseBiologyModule {
      rodCell() : BaseBiologyModule(gAllBmEvents){}

      template <typename T>
      void Run(T* cell){
        //what will make the cell grow
        if (cell->GetCellType() > 0) {
          array<double, 3> direction {0, 1, 0};
          //make it grow along the y axis
          cell->UpdateMassLocation(direction);
          //to ensure that there is no traction between cells
          cell->SetTractorForce({0,0,0});
          cell->SetCellType(cell->GetCellType()-1);
        }
      }
    };

//may not need this at all --> just extra function to define the substance?
  // struct GaussianBand {
  //   double mean_;
  //   double sigma_;
  //   uint8_t axis_;
  //
  //   GaussianBand(double mean, double sigma, uint8_t axis) {
  //     mean_ = mean;
  //     sigma_ = sigma;
  //     axis_ = axis;
  //   }
  //
  //   double operator()(double x, double y, double z) {
  //     switch(axis_) {
  //       case Axis::kXAxis: return ROOT::Math::normal_pdf(x, sigma_, mean_);
  //       case Axis::kYAxis: return ROOT::Math::normal_pdf(y, sigma_, mean_);
  //       case Axis::kZAxis: return ROOT::Math::normal_pdf(z, sigma_, mean_);
  //       default: throw std::logic_error("You have chosen an non-existing axis!");
  //     }
  //   }
  // };

// Define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  /*using BiologyModules = Variant<Chemotaxis, kSecretion, midgetCell,
  parasolCell, bistratifiedCell, amacrineCell, bipolarCell,
  horizontalCell, coneCell, rodCell>;*/
  using BiologyModules = Variant<midgetCell,
  parasolCell, bistratifiedCell, amacrineCell, bipolarCell,
  horizontalCell, coneCell, rodCell, Chemotaxis>;
  using AtomicTypes = VariadicTypedef <MyCell>;
};

inline int Simulate(int argc, const char** argv) {
  InitializeBioDynamo(argc, argv);

  /* this is to ensure that the model is reproducible for a
  specific seed
  part 2: TODO check what is a good value range for a seed to occur
  */
    int randSeed = rand() % 1000;
    gTRandom.SetSeed(randSeed);
    cout << "modelling seed: " << randSeed <<endl;

  /* Define initial model - in this example: single cell at origin
  <s>Assume if we start from outer to inner retina --> start with RPE</s>

  Part 2: should instead work from ganglion cells downwards
  assume how the light rays fall onto the retina
  */

  /*rods
  part 2:
  cells are longish -> outer seg + inner seg + nucleus
  connect to horizontal cells + bipolar cells + Parasol cells? <check this>
  */
    auto construct_rod = [](const std::array<double, 3>& position){
      MyCell cell(position);
      cell.SetDiameter(20);
      cell.AddBiologyModule(rodCell());
      cell.AddBiologyModule(Chemotaxis());
      cell.SetCellType(8);
      return cell;
    };
    cout << "Rod cells created" << endl;

    //this creates the cell at random: min bound, max bound, no of cell, cell type
     ModelInitializer::CreateCellsRandom(350, gTRandom.Uniform(350,Param::max_bound_),
                                       20, construct_rod);

  /*cones
  part 2:
  cells are long but rounder and wider than rods
  connect to horizontal cells + bipolar cells + Parasol cells? <check this>
  */

    auto construct_cone = [](const std::array<double, 3>& position){
      MyCell cell(position);
      cell.SetDiameter(20);
      cell.AddBiologyModule(coneCell());
      cell.AddBiologyModule(Chemotaxis());
      cell.SetCellType(7);
      return cell;
    };
    cout << "Cone cells created" << endl;

    //this creates the cell at random: min bound, max bound, no of cell, cell type
     ModelInitializer::CreateCellsRandom(350, gTRandom.Uniform(350,Param::max_bound_),
                                       20, construct_cone);
/*horizontal cells
part 2:
cells work laterally
connected to outputs from rods and cones
*/

    auto construct_horizontal = [](const std::array<double, 3>& position){
      MyCell cell(position);
      cell.SetDiameter(30);
      cell.AddBiologyModule(horizontalCell());
      cell.AddBiologyModule(Chemotaxis());
      cell.SetCellType(6);
      return cell;
    };
    cout << "Horizontal cells created" << endl;

    //this creates the cell at random: min bound, max bound, no of cell, cell type
     ModelInitializer::CreateCellsRandom(300, gTRandom.Uniform(300,450),
                                       20, construct_horizontal);

  /*bipolar cells
  part 2: exists between photoreceptors and ganglion cells
  so it will either synpase with:
  1. photoreceptors -> rods/ cones (i.e. Parasol cells)
  2. horizontal cells
  */
    auto construct_bipolar = [](const std::array<double, 3>& position){
      MyCell cell(position);
      cell.SetDiameter(40);
      cell.AddBiologyModule(bipolarCell());
      cell.AddBiologyModule(Chemotaxis());
      cell.SetCellType(5);
      return cell;
    };
    cout << "Bipolar cells created" << endl;

    //this creates the cell at random: min bound, max bound, no of cell, cell type
     ModelInitializer::CreateCellsRandom(200, gTRandom.Uniform(200,400),
                                       20, construct_bipolar);

  /*amacrine cells
  part 2: This is in the inner limiting layer where it links to the
  bipolar and ganglion cells
  so it takes the input from the ganglion cell o the bipolar cell
  So for ganglion cells --> link to the MIDGET & BISTRATIFIED CELLS
  Cells work laterally
  */
    auto construct_amacrine = [](const std::array<double, 3>& position){
      MyCell cell(position);
      cell.SetDiameter(45);
      cell.AddBiologyModule(amacrineCell());
      cell.AddBiologyModule(Chemotaxis());
      cell.SetCellType(4);
      return cell;
    };
    cout << "Amacrine cells created" << endl;

    //this creates the cell at random: min bound, max bound, no of cell, cell type
     ModelInitializer::CreateCellsRandom(100, gTRandom.Uniform(100,400),
                                       20, construct_amacrine);

/* ganglion cell
part 2: there are various types of ganglion cells to take into account that
connects to different cells for inputs
for now, take into consideration the following:
1. Midget cell -> input from BIPOLAR CELLS + small size
2. Parasol cell -> input from rod & cones; part of M pathway + large size (dentritic)
3. Bistratified cell -> input from BIPOLAR & AMACRINE CELLS + medium size
4. Photosensitive ganglion cells <ignore>
5. other cells for saccades <ignore>
create a condition to link to the layer below it -> "if statements"?
*/

    auto construct_midget = [](const std::array<double, 3>& position){
      MyCell cell(position);
      cell.SetDiameter(50);
      cell.AddBiologyModule(midgetCell());
      cell.AddBiologyModule(Chemotaxis());
      cell.SetCellType(3);
      return cell;
    };
    cout << "Midget cells created" << endl;

    //this creates the cell at random: min bound, max bound, no of cell, cell type
     ModelInitializer::CreateCellsRandom(Param::min_bound_, gTRandom.Uniform(0,200),
                                       20, construct_midget);
    //defining substances in simulation
    //diffusion coefficient of 0.5, a decay constant 0f 0.1 and a resolution of 1
  //  ModelInitializer::DefineSubstance(kSubstance, "Substance", 0.5, 0.1, 1);

    auto construct_parasol = [](const std::array<double, 3>& position){
      MyCell cell(position);
      cell.SetDiameter(100);
      cell.AddBiologyModule(parasolCell());
      cell.AddBiologyModule(Chemotaxis());
      cell.SetCellType(2);
      return cell;
    };
    cout << "Parasol cells created" << endl;

    //this creates the cell at random: min bound, max bound, no of cell, cell type
     ModelInitializer::CreateCellsRandom(Param::min_bound_, Param::max_bound_,
                                       20, construct_parasol);

    // for (int i=0; i<20; i++) {
    //   array<double, 3> position = {gTRandom.Uniform(0,500), gTRandom.Uniform(0,500), 0};
    //   MyCell cell(position);
    //   cell.SetDiameter(75);
    //   cell.AddBiologyModule(bistratifiedCell());
    //   cell.AddBiologyModule(Chemotaxis());
    //   cell.SetCellType(90);
    //   ResourceManager<>::Get()->push_back(cell);
    //   if (position[0] == 0 && position[1] == 0 && position[2] == 0) {
    //     cell.AddBiologyModule(kSecretion());
    //   }
    // }

    auto construct_bistratified = [](const std::array<double, 3>& position) {
      MyCell cell(position);
      cell.SetDiameter(75);
      cell.AddBiologyModule(bistratifiedCell());
      cell.AddBiologyModule(Chemotaxis());
      cell.SetCellType(1);
      return cell;
    };
    //CellCreator(Param::min_bound_, Param::max_bound_, 50, construct_bistratified);
    cout << "Bistratified cells created" << endl;

   //this creates the cell at random: min bound, max bound, no of cell, cell type
   //TODO change min and max bound so that cells are created within a
   //certain layer --> should i do gTRandom.Uniform(0,200) ???
    ModelInitializer::CreateCellsRandom(Param::min_bound_, gTRandom.Uniform(0,200),
                                      20, construct_bistratified);
    //defining substances in simulation
    //diffusion coefficient of 0.5, a decay constant 0f 0.1 and a resolution of 1
    ModelInitializer::DefineSubstance(kSubstance, "Substance", 0.5, 0.1, 1);
    //initialise substance: enum of substance, name, function type used
    //mean value of 120 along the x-axis, and a variance of 5
    //along which axis (0 = x, 1 = y, 2 = z). See the documentation of `GaussianBand` for
    //information about its arguments
    ModelInitializer::InitializeSubstance(kSubstance, "Substance",
      GaussianBand(120, 5, Axis::kXAxis));


  //link to paraview to show visualization
    Param::live_visualization_ = true;
    Param::min_bound_ = 0;
    Param::max_bound_ = 520;
    Param::run_mechanical_interactions_ = true;

  // Run simulation for one timestep
  Scheduler<> scheduler;
  int maxStep = 1000;
  for (int i = 0; i < maxStep; i++){
    scheduler.Simulate(1);
    //for every 10 cells
    // if(i %10==0){
    //
    // }
  }


  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

} // namespace bdm

#endif // RETINA_PROJ2_H_
