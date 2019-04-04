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
          rreg->min(d) = Eval::evaluate( head->subHeads[d]->params[0]->partA, variables );
          rreg->max(d) = Eval::evaluate( head->subHeads[d]->params[1]->partA, variables );
        }
        region = rreg;
      }
    }

    return region;
  } 

  void ParseConstructor::parse_particle_template(HeadNode *head, const std::map<string, string>& variables, std::map<string, ParticleTemplate>& particle_templates, GFlow *gflow) {
    // Create a particle template to set
    ParticleTemplate p_template;

    ParseHelper parser(head);
    parser.set_variables(variables);
    // Declare valid options
    parser.addValidSubheading("Radius");
    parser.addValidSubheading("Mass");
    parser.addValidSubheading("Type");
    // Make sure only valid options were used
    if (!parser.checkValidSubHeads()) {
      cout << "Warning: Invalid Headings:\n";
      for (auto ih : parser.getInvalidSubHeads())
        cout << " -- Heading: " << ih << endl;
    }
    // Sort options
    parser.sortOptions();
    // Pointer for head nodes
    HeadNode *hd = nullptr;
    
    // --- Look for options
    string option, type;
    
    // --- Look for Radius option
    parser.getHeading_Necessary("Radius");
    hd = parser.first();
    p_template.radius_engine = getRandomEngine(hd, variables, type, gflow);
    p_template.radius_string = type;

    // --- Look for Mass option
    parser.getHeading_Necessary("Mass");
    hd = parser.first();
    parser.extract_first_parameter(option, hd);
    RealType m;
    if (option=="Density" && parser.extract_first_arg(m)) {
      p_template.mass_engine = new DeterministicEngine(m);
      p_template.mass_string = "Density";
    }
    else if (parser.extract_first_arg(m)) {
      p_template.mass_engine = new DeterministicEngine(m);
      p_template.mass_string = "Mass";
    }

    // --- Look for Type option
    parser.getHeading_Necessary("Type");
    hd = parser.first();
    p_template.type_engine = getRandomEngine(hd, variables, type, gflow);
    p_template.type_string = type;

    // Add to particle templates
    particle_templates.insert(std::pair<string, ParticleTemplate>(head->params[0]->partA, p_template));
  }

  RandomEngine* ParseConstructor::getRandomEngine(HeadNode *h, const std::map<string, string>& variables, string &type, GFlow *gflow) {
    // Clear type string
    type.clear();

    ParseHelper parser(h);
    parser.set_variables(variables);

    // Check structure
    if (h->params.size()!=1)
      throw BadStructure("Random engine needs one parameter, found "+toStr(h->params.size())+".");
    else if (!h->params[0]->partB.empty()) {
      // A type is specified
      type = h->params[0]->partA;
      // We expect a number in partB
      return new DeterministicEngine(parser.value<RealType>(h->params[0]->partB));
    }
    // If there is no body
    else if (h->subHeads.empty()) {
      string token = h->params[0]->partA;
      // We expect either a number in partA, or a string (e.g. "Inf").
      if (isdigit(token.at(0))) {
        return new DeterministicEngine(parser.value<RealType>(token));
      }
      else if (h->heading=="Type" && token=="Equiprobable") {
        type = token;
        int NTypes = gflow->getNTypes();
        // Only for type - Create each type with equal probability
        RealType p = 1./NTypes;
        vector<RealType> probabilities, values;
        for (int i=0; i<NTypes; ++i) {
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
    string param = h->params[0]->partA;
    // Check for options:
    if (param=="Random") {
      // Discrete Random, with specified values and probabilities
      vector<RealType> probabilities, values;
      for (auto sh : h->subHeads) {
        if (sh->heading!="P")
          throw BadStructure("Discrete probability should be indicated with a 'P.'");
        if (sh->params.size()!=2)
          throw BadStructure("Expect value and probability, encountered "+toStr(sh->params.size())+" options.");
        if (!sh->params[0]->partB.empty() || !sh->params[1]->partB.empty())
          throw BadStructure("Expected one part parameters in discrete probability specification.");
        // Store the value and its probability
        values.push_back(parser.value<RealType>(sh->params[0]->partA)); 
        probabilities.push_back(parser.value<RealType>(sh->params[1]->partA));
      }
      return new DiscreteRandomEngine(probabilities, values);
    }
    else if (param=="Uniform") {
      if (h->subHeads.size()!=2) 
        throw BadStructure("Uniform distribution needs a min and a max. Found "+toStr(h->subHeads.size())+" parameters.");
      HeadNode *lwr = h->subHeads[0], *upr = h->subHeads[1];
      // Swap if necessary
      if (lwr->heading!="Min") std::swap(lwr, upr);
      // Check structure
      if (lwr->heading!="Min" || upr->heading!="Max") 
        throw BadStructure("Headings incorrect for uniform distribution. Need [Min] and [Max]");
      if (lwr->params.size()!=1 || upr->params.size()!=1) 
        throw BadStructure("More than one parameter for uniform distribution.");
      // Set the uniform random engine
      return new UniformRandomEngine(
        parser.value<RealType>(lwr->params[0]->partA), parser.value<RealType>(upr->params[0]->partA)
      );
    }
    else if (param=="Normal") {
      if (h->subHeads.size()!=2) 
        throw BadStructure("Normal distribution needs average and variance. Found "+toStr(h->subHeads.size())+" parameters.");
      HeadNode *ave = h->subHeads[0], *var = h->subHeads[1];
      // Swap if necessary
      if (ave->heading!="Ave") std::swap(ave, var);
      // Check structure
      if (ave->heading!="Ave" || var->heading!="Var")
        throw BadStructure("Headings incorrect for normal distribution.");
      if (ave->params.size()!=1 || var->params.size()!=1) 
        throw BadStructure("More than one parameter for normal distribution.");
      // Set the uniform random engine
      return new NormalRandomEngine(
        parser.value<RealType>(ave->params[0]->partA), parser.value<RealType>(var->params[0]->partA)
      );
    }
    else throw BadStructure("Unrecognized choice for a random engine.");
    // We should never reach here
    throw BadStructure("An unreachable part of code was reached!");
    // Token return
    return nullptr;
  }

}