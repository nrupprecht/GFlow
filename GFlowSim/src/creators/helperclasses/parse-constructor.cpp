#include "parse-constructor.hpp"

namespace GFlowSimulation {

  Region* ParseConstructor::parse_region(HeadNode *head, const std::map<string, string> & variables, GFlow *gflow) {
    // Pointer to a region
    Region *region = nullptr;
    // Simulation dimensions
    int sim_dimensions = gflow->getSimDimensions();

    // Check if we should use the full simulation bounds - the bounds are then rectangular
    if (head->subHeads.empty() && head->params.size()>0 && head->params[0]->partA=="Full")
      region = new RectangularRegion(sim_dimensions, gflow->getBounds());
    // Spherical bounds
    else if (head->params.size()>0 && head->params[0]->partA=="Sphere") {
      if (head->subHeads.size()!=2)
        throw BadStructure("Expected position and radius for spherical bounds.");
      if (head->subHeads[0]->params.size()!=sim_dimensions)
        throw BadStructure("Vector needs the correct number of dimensions.");
      SphericalRegion *sreg = new SphericalRegion(sim_dimensions);
      // Get the center
      for (int d=0; d<sim_dimensions; ++d)
        sreg->center(d) = Eval::evaluate(head->subHeads[0]->params[d]->partA, variables);
      // Get the radius
      sreg->radius() = Eval::evaluate(head->subHeads[1]->params[0]->partA, variables);
      region = sreg;
    }
    // Otherwise, rectangular bounds
    else {
      // If the bounds is a Bounds: Box=[width]
      if (head->params.size()>0 && head->params[0]->partA=="Box") {
        if (head->subHeads.size()!=0)
          throw BadStructure("Box bounds should not have any subheads");
        if (head->params.size()!=1)
          throw BadStructure("Box bounds needs one parameter, found "+toStr(head->params.size())+".");
        RectangularRegion *rreg = new RectangularRegion(sim_dimensions);
        RealType width = Eval::evaluate( head->params[0]->partB, variables );
        // Set all min/max to be 0/width
        for (int d=0; d<sim_dimensions; ++d) {
          rreg->min(d) = 0;
          rreg->max(d) = width;
        }
        region = rreg;
      }
      // Otherwise, all bounds min/max are specified.
      else {
        if (head->subHeads.size()!=sim_dimensions) 
          throw BadStructure("Expected "+toStr(sim_dimensions)+" arguments, found "+toStr(head->subHeads.size()));
        // Set bounds
        RectangularRegion *rreg = new RectangularRegion(sim_dimensions);
        for (int d=0; d<head->subHeads.size(); ++d) {
          // Check for well formed options.
          if (head->subHeads[d]->params.size()!=2 || !head->subHeads[d]->subHeads.empty()) 
            throw BadStructure("Bounds need a min and a max, we found "+toStr(head->subHeads[d]->params.size())+" parameters.");
          // Extract the bounds.
	  RealType mn = Eval::evaluate( head->subHeads[d]->params[0]->partA, variables );
	  RealType mx = Eval::evaluate( head->subHeads[d]->params[1]->partA, variables );
          rreg->min(d) = mn;
          rreg->max(d) = mx;
        }
        region = rreg;
      }
    }

    return region;
  } 

  void ParseConstructor::parse_particle_template(HeadNode *head, const std::map<string, string>& variables, std::map<string, ParticleTemplate>& particle_templates, GFlow *gflow) {
    // Create a particle template to set
    ParticleTemplate p_template;

    // Create parser.
    TreeParser parser(head, variables);
    parser.addHeadingNecessary("Radius", "We need radius info for particle templates!");
    parser.addHeadingNecessary("Mass", "We need mass info for particle templates!");
    parser.addHeadingNecessary("Type", "We need type info for particle templates!");
    // Check headings for validity.
    parser.check();

    // --- Look for options
    string option, type;

    // --- Look for Radius option
    parser.focus0("Radius");
    p_template.radius_engine = getRandomEngine(parser.getNode(), variables, type, gflow);
    p_template.radius_string = type;
    
    // --- Look for Mass option
    parser.focus0("Mass");
    RealType m = 0;
    // Density
    if (parser.argName()=="Density") {
      // Extract the density value
      if (!parser.val(m)) throw BadStructure("Mass specified as DENSITY, but value could not be extracted.");
      // Set the mass engine
      p_template.mass_engine = new DeterministicEngine(m);
      p_template.mass_string = "Density";
    }
    // Mass
    else if (parser.val(m)) {
      p_template.mass_engine = new DeterministicEngine(m);
      p_template.mass_string = "Mass";
    }
    // Error
    else throw BadStructure("Mass specification not recognized.");

    // --- Look for Type option
    parser.focus0("Type");
    p_template.type_engine = getRandomEngine(parser.getNode(), variables, type, gflow);
    p_template.type_string = type;

    // Add to particle templates
    particle_templates.insert(std::pair<string, ParticleTemplate>(head->params[0]->partA, p_template));
  }

  RandomEngine* ParseConstructor::getRandomEngine(HeadNode *h, const std::map<string, string>& variables, string &type, GFlow *gflow) {
    // Clear type string
    type.clear();

    // Create a parser
    TreeParser parser(h, variables);

    // Check structure
    if (parser.args_size()!=1)
      throw BadStructure("Random engine needs one parameter, found "+toStr(h->params.size())+".");
    else if (parser.valName()!="") {
      // A type is specified
      type = parser.argName();
      // We expect a number as the value
      RealType v;
      if (parser.arg(v)) return new DeterministicEngine(v);
      else throw BadStructure("Expected a value for Deterministic Engine.");
    }
    // If there is no body
    else if (parser.body_size()==0) {
      // This will allow us to recognize variables.
      string token;
      parser.arg(token);

      // We expect either a number in partA, or a string (e.g. "inf").
      if (isdigit(token.at(0))) {
        RealType v;
        if (parser.arg(v)) return new DeterministicEngine(v);
        else throw BadStructure("Expected a value for DeterministicEngine.");
      }
      else if (parser.heading()=="Type" && token=="Equiprobable") {
        type = token;
        // Only for type - Create each type with equal probability
        int ntypes = gflow->getNTypes();
        RealType p = 1./ntypes;
        vector<RealType> probabilities, values;
        for (int i=0; i<ntypes; ++i) {
          probabilities.push_back(p);
          values.push_back(i);
        }
        return new DiscreteRandomEngine(probabilities, values);
      }
      else {
        type = token;
        // The engine will not matter
        return new DeterministicEngine(0);
      }
    }

    // Parameter string
    string param = parser.argName();
    // Check for options:
    if (param=="Random") {
      // Discrete Random, with specified values and probabilities
      vector<RealType> probabilities, values;

      if (parser.begin()) {
        do {
          if (parser.heading()!="P") 
             throw BadStructure("Discrete probability should be indicated with a 'P.'");

          if (parser.args_size()!=2)
            throw BadStructure("Expect value and probability, encountered "+toStr(parser.args_size())+" options.");
          if (parser.valName(0)!="" || parser.valName(1)!="")
            throw BadStructure("Expected one part parameters in discrete probability specification.");
          // Store the value and its probability
          RealType v = 0, p = 0;
          // If the value, probability are there, store them
          if (parser.arg(v, 0) && parser.arg(p, 0)) {
            values.push_back(v); 
            probabilities.push_back(p);
          }
          // Otherwise, throw an exception
          else throw BadStructure("Expected value, probability. One or both were missing.");
        } while (parser.next());
        // Return the random engine
        return new DiscreteRandomEngine(probabilities, values);
      }
      else throw BadStructure("Parser couldn't set up iteration.");
    }
    else if (param=="Uniform") {
      if (parser.body_size()!=2) 
        throw BadStructure("Uniform distribution needs a min and a max. Found "+toStr(parser.body_size())+" parameters.");

      RealType mn = 0, mx = 0;
      // Set the uniform random engine if we have the data we need.
      if (parser.firstArg("Min", mn) && parser.firstArg("Max", mx)) return new UniformRandomEngine(mn, mx);
      // Otherwise, throw exception.
      else throw BadStructure("Needed Min, Max, couldn't find one or both of those.");
    }
    else if (param=="Normal") {
      if (parser.body_size()!=2) 
        throw BadStructure("Normal distribution needs 2 parameters. Found "+toStr(parser.body_size())+" parameters.");
      // Extract average, variance
      RealType ave = 0, var = 0;
      // Set up random engine
      if (parser.firstArg("Ave", ave) && parser.firstArg("Var", var)) return new NormalRandomEngine(ave, var);
      else throw BadStructure("Normal distribution needs Ave, Var, couldn't find one or both of those.");
    }
    else throw BadStructure("Unrecognized choice for a random engine.");
    // We should never reach here
    throw BadStructure("An unreachable part of code was reached!");
    // Token return
    return nullptr;
  }

}
