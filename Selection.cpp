#include "rts/operator/Selection.hpp"
#include "rts/operator/PlanPrinter.hpp"
#include "rts/database/Database.hpp"
#include "rts/runtime/Runtime.hpp"
#include "rts/segment/DictionarySegment.hpp"
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#ifdef __GNUC__
#if (__GNUC__>4)||((__GNUC__==4)&&(__GNUC_MINOR__>=5))
#define CONFIG_TR1
#endif
#endif
#ifdef CONFIG_TR1
#include <tr1/regex>
#endif
//---------------------------------------------------------------------------
// RDF-3X
// (c) 2008 Thomas Neumann. Web site: http://www.mpi-inf.mpg.de/~neumann/rdf3x
//
// This work is licensed under the Creative Commons
// Attribution-Noncommercial-Share Alike 3.0 Unported License. To view a copy
// of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/
// or send a letter to Creative Commons, 171 Second Street, Suite 300,
// San Francisco, California, 94105, USA.
//---------------------------------------------------------------------------
using namespace std;
//---------------------------------------------------------------------------
Selection::RetrieveGeoIDsIndexScan::Hint::Hint(RetrieveGeoIDsIndexScan& scan)
   : scan(scan)
   // Constructor
{
}
//---------------------------------------------------------------------------
Selection::RetrieveGeoIDsIndexScan::Hint::~Hint()
   // Destructor
{
}
//---------------------------------------------------------------------------
void Selection::RetrieveGeoIDsIndexScan::Hint::next(unsigned& value1,unsigned& value2,unsigned& value3)
	// Do Nothing
{
}
//---------------------------------------------------------------------------
Selection::RetrieveGeoIDsIndexScan::RetrieveGeoIDsIndexScan(Database& db, Database::DataOrder order)
   : facts(db.getFacts(order)), order(order), scan(0), hint(*this)
   // Constructor
{
}
//---------------------------------------------------------------------------
Selection::RetrieveGeoIDsIndexScan::~RetrieveGeoIDsIndexScan()
   // Destructor
{
}
//---------------------------------------------------------------------------
unsigned Selection::RetrieveGeoIDsIndexScan::getGeoID(unsigned filter1, unsigned filter2) 
{
	if (!this->scan.first(this->facts,filter1,filter2,0)) {
		cout << "SCAN FIRST ERROR!!" << endl;	
		exit(-1);
	}

	if ((this->scan.getValue1()>filter1)||((this->scan.getValue1()==filter1)&&(this->scan.getValue2()>filter2))) {
		cout << "SCAN GET ENTRY ERROR!!" << endl;
		exit(-1);
	}	

	return (this->scan.getValue3());
}
//---------------------------------------------------------------------------
Selection::Result::~Result()
   // Destructor
{
}
//---------------------------------------------------------------------------
bool Selection::DistanceF::getCellId(unsigned value, unsigned &c_id, unsigned &g_pos, unsigned &level)
	// Extract the cell id, the level, and the grid pos from the ID of an entity. 
	// Return true if the input entity is associated with spatial info, false otherwise.
{
	//cout << "Given value: " << value << endl;

	//XXX: attention attention
	//unsigned mask = 30;
	//unsigned mask = 15;

	level = 0;
	c_id = value;

	/*
	while((c_id & (unsigned) 1) && (level<MAX_LEVEL)){

		c_id >>= 1;
		level++;		
	}

	if(level)
		c_id >>= ((b1+b2)-2*level-level);
	*/

	if(c_id & 1){
		level = ((c_id & LVMASK) >> 1);
	
		//XXX: attention attention
		level += LVPREF;
	}
	else return false;

	c_id >>= (grid->b-2*level);

	//cout << "Cell Id: " << c_id << endl;
	//cout << "level: " << level << endl;
	//getchar();

	if(!level) 
		c_id = g_pos = (grid->size-1);
	else{
		g_pos = c_id;
		unsigned step = grid->MAX_LEVEL;
		while(step>level){
		
			g_pos += ((unsigned) 1<<(grid->b1-2*(grid->MAX_LEVEL-step)));		//pow(2,b1-2*(MAX_LEVEL-step));
			step--;
		} 
	}

	//cout << "Cell Id: " << c_id << endl;
	//cout << "g_pos: " << g_pos << endl;
	//getchar();

	//XXX: conflicts with discretized coordinates must be checked here
	//if(grid->grid[g_pos]>value){
	//if(value & ((unsigned) 1)){//if the LSB is set

		//cout << "Spatial." << endl;
		return true;
	//}

	//the entity is not associated with spatial information
	//cout << "Not spatial." << endl;
	//return false;

}
//---------------------------------------------------------------------------
bool Selection::WithinF::getCellId(unsigned value, unsigned &c_id, unsigned &g_pos, unsigned &level)
	// Extract the cell id, the level, and the grid pos from the ID of an entity. 
	// Return true if the input entity is associated with spatial info, false otherwise.
{
	//cout << "Given value: " << value << endl;

	//XXX: attention attention
	//unsigned mask = 30;
	//unsigned mask = 15;

	level = 0;
	c_id = value;

	/*
	while((c_id & (unsigned) 1) && (level<MAX_LEVEL)){

		c_id >>= 1;
		level++;		
	}

	if(level)
		c_id >>= ((b1+b2)-2*level-level);
	*/

	if (c_id & 1) { 
		level = ((c_id & LVMASK) >> 1);
		
		//XXX: attention attention
		#if GRID
		if (level >= abs(LVPREF))
		#endif
		level += LVPREF;
	}
	else return false;
	
	c_id >>= (grid->b - 2*level);

	//cout << "Cell Id: " << c_id << endl;
	//cout << "level: " << level << endl;
	//getchar();

	if (!level) 
		c_id = g_pos = (grid->size - 1);
	else {
		g_pos = c_id;
		unsigned step = grid->MAX_LEVEL;
		while (step > level) {
			g_pos += ((unsigned) 1 << (grid->b1 - 2*(grid->MAX_LEVEL - step)));  // pow(2,b1-2*(MAX_LEVEL-step));
			step--;
		} 
	}

	//cout << "Cell Id: " << c_id << endl;
	//cout << "g_pos: " << g_pos << endl;
	//getchar();

	//XXX: conflicts with discretized coordinates must be checked here
	//if(grid->grid[g_pos]>value){
	//if(value & ((unsigned) 1)){//if the LSB is set

		//cout << "Spatial." << endl;
		return true;
	//}

	//the entity is not associated with spatial information
	//cout << "Not spatial." << endl;
	//return false;
}
//---------------------------------------------------------------------------
bool Selection::DistanceF::evalPredicate(double minDist, double maxDist, bool& verified)
{
	switch (type)
	{
		case QueryGraph::Filter::Less: 
		{
			if(minDist>=distThreshold){
				//std::cout << "DistanceF: Filtered." << std::endl;
				return false;
			}
			else if(maxDist<distThreshold){	
				verified = true;
				//std::cout << "DistanceF: Verified." << std::endl;
				return true;
			}
			
			verified = false;
			return true;
		}
		case QueryGraph::Filter::LessOrEqual: 
		{
			if(minDist>distThreshold)
				return false;
			else if(maxDist<=distThreshold){	
				verified = true;
				return true;
			}
			
			verified = false;
			return true;
		}
		case QueryGraph::Filter::Greater: 
		{
			if(minDist>distThreshold){
				verified = true;
				return true;
			}
			else if(maxDist<=distThreshold)
				return false;
			
			verified = false;
			return true;
		}
		case QueryGraph::Filter::GreaterOrEqual:
		{
			if(minDist>=distThreshold){
				verified = true;
				return true;
			}
			else if(maxDist<distThreshold)
				return false;
			
			verified = false;
			return true;	
		}
		case QueryGraph::Filter::Equal:
		{
			if(minDist>distThreshold || maxDist<distThreshold)
				return false;

			verified = false;
			return true;				
		}  
		default:
		{
			std::cout << "Unknown spatial join predicate!!" << std::endl;		
			return true;
		}
	}
}
//---------------------------------------------------------------------------
void Selection::DistanceF::getCenterCoords(unsigned cord_x1, unsigned cord_y1, unsigned cord_x2, unsigned cord_y2, double *c)
{
	if (cord_x1 == cord_x2) {  // cell refers to the max level
		c[0] = cord_x1 + 0.5;
		c[1] = cord_y1 + 0.5;
	}
	else {  // cell refers to a higher level
		c[0] = cord_x1 + ((double) (cord_x2 + 1 - cord_x1)) / 2;
		c[1] = cord_y1 + ((double) (cord_y2 + 1 - cord_y1)) / 2;
	}		
}
//---------------------------------------------------------------------------
void Selection::DistanceF::getMinMaxDist(unsigned cord_x1, unsigned cord_y1, unsigned cord_x2, unsigned cord_y2, unsigned c_x1, 
													  unsigned c_y1, unsigned c_x2, unsigned c_y2, double &minDist, double &maxDist)
{
	double c1[2], c2[2], coord1[2], coord2[2];
	// get center coordinates
	getCenterCoords(cord_x1, cord_y1, cord_x2, cord_y2, c1);
	getCenterCoords(c_x1, c_y1, c_x2, c_y2, c2);

	if (c1[0] <= c2[0]) {  // x<=x'

		if (c1[1] <= c2[1]) {  // y<=y'

			// left cell is south west of the right one

			bool f1 = ((cord_x2 + 1) >= c_x1), f2 = ((cord_y2 + 1) >= c_y1);

			if (f1 && f2) {  // they meet at least
				minDist = 0;
			}
			else if (f1) {  // distance in y (x must be common)
				getUpperLeftPoint(cord_x2, cord_y2, coord1);
				getLowerLeftPoint(cord_x2, c_y1, coord2);
				//coord1[0] = cord_x2; coord1[1] = cord_y2+1;
				//coord2[0] = cord_x2; coord2[1] = c_y1;
				minDist = dist2(coord1, coord2);
			}
			else if (f2) {  // distance in x (y must be common)
				getLowerRightPoint(cord_x2, c_y1, coord1);
				getLowerLeftPoint(c_x1, c_y1, coord2);
				//coord1[0] = cord_x2+1; coord1[1] = c_y1;
				//coord2[0] = c_x1; coord2[1] = c_y1;
				minDist = dist2(coord1, coord2);
			}
			else { 
				getUpperRightPoint(cord_x2, cord_y2, coord1);
				getLowerLeftPoint(c_x1, c_y1, coord2);
				//coord1[0] = cord_x2+1; coord1[1] = cord_y2+1;
				//coord2[0] = c_x1; coord2[1] = c_y1;
				minDist = dist2(coord1, coord2);				
			}

			getLowerLeftPoint(cord_x1, cord_y1, coord1);
			getUpperRightPoint(c_x2, c_y2, coord2);
			//coord1[0] = cord_x1; coord1[1] = cord_y1;
			//coord2[0] = c_x2+1; coord2[1] = c_y2+1;
			maxDist = dist2(coord1, coord2);

			return;
		}
		else {  // y>y'

			// left cell is north west of the right one

			bool f1 = ((cord_x2 + 1) >= c_x1), f2 = (cord_y1 <= (c_y2 + 1));

			if (f1 && f2) {  // they meet at least
				minDist = 0;
			}
			else if (f1) {  // distance in y (x must be common)
				getLowerLeftPoint(cord_x2, cord_y1, coord1);
				getUpperLeftPoint(cord_x2, c_y2, coord2);
				//coord1[0] = cord_x2; coord1[1] = cord_y1;
				//coord2[0] = cord_x2; coord2[1] = c_y2+1;
				minDist = dist2(coord1, coord2);
			}
			else if (f2) {  // distance in x (y must be common)
				getLowerRightPoint(cord_x2, cord_y2, coord1);
				getLowerLeftPoint(c_x1, cord_y2, coord2);
				//coord1[0] = cord_x2+1; coord1[1] = cord_y2;
				//coord2[0] = c_x1; coord2[1] = cord_y2;
				minDist = dist2(coord1, coord2);
			}
			else {
				getLowerRightPoint(cord_x2, cord_y1, coord1);
				getUpperLeftPoint(c_x1, c_y2, coord2);
				//coord1[0] = cord_x2+1; coord1[1] = cord_y1;
				//coord2[0] = c_x1; coord2[1] = c_y2+1;
				minDist = dist2(coord1, coord2);				
			}

			getUpperLeftPoint(cord_x1, cord_y2, coord1);
			getLowerRightPoint(c_x2, c_y1, coord2);
			//coord1[0] = cord_x1; coord1[1] = cord_y2+1;
			//coord2[0] = c_x2+1; coord2[1] = c_y1;	
			maxDist = dist2(coord1, coord2);

			return;
		}

	}
	else if (c1[1] <= c2[1]) {  // x>x' && y<=y'

		// left cell is south east of the right one

		bool f1 = (cord_x1 <= (c_x2 + 1)), f2 = ((cord_y2 + 1) >= c_y1);

		if (f1 && f2) {  // they meet at least
			minDist = 0;
		}
		else if (f1) {  // distance in y (x must be common)
			getUpperLeftPoint(cord_x1, cord_y2, coord1);
			getLowerLeftPoint(cord_x1, c_y1, coord2);
			//coord1[0] = cord_x1; coord1[1] = cord_y2+1;
			//coord2[0] = cord_x1; coord2[1] = c_y1;
			minDist = dist2(coord1, coord2);
		}
		else if (f2) {  // distance in x (y must be common)
			getLowerRightPoint(c_x2, c_y1, coord1);
			getLowerLeftPoint(cord_x1, c_y1, coord2);
			//coord1[0] = c_x2+1; coord1[1] = c_y1;
			//coord2[0] = cord_x1; coord2[1] = c_y1;
			minDist = dist2(coord1, coord2);
		}
		else {
			getLowerRightPoint(c_x2, c_y1, coord1);
			getUpperLeftPoint(cord_x1, cord_y2, coord2);
			//coord1[0] = c_x2+1; coord1[1] = c_y1;
			//coord2[0] = cord_x1; coord2[1] = cord_y2+1;
			minDist = dist2(coord1, coord2);				
		}

		getUpperLeftPoint(c_x1, c_y2, coord1);
		getLowerRightPoint(cord_x2, cord_y1, coord2);
		//coord1[0] = c_x1; coord1[1] = c_y2+1;
		//coord2[0] = cord_x2+1; coord2[1] = cord_y1;
		maxDist = dist2(coord1, coord2);

		return;
	}
	else {  // x>x' && y>y'

		// left cell is north east of the right one

		bool f1 = (cord_x1 <= (c_x2 + 1)), f2 = (cord_y1 <= (c_y2 + 1));

		if (f1 && f2) {  // they meet at least
			minDist = 0;
		}
		else if (f1) {  // distance in y (x must be common)
			getLowerLeftPoint(cord_x1, cord_y1, coord1);
			getUpperLeftPoint(cord_x1, c_y2, coord2);
			//coord1[0] = cord_x1; coord1[1] = cord_y1;
			//coord2[0] = cord_x1; coord2[1] = c_y2+1;
			minDist = dist2(coord1, coord2);
		}
		else if (f2) {  // distance in x (y must be common)
			getLowerLeftPoint(cord_x1, cord_y1, coord1);
			getLowerRightPoint(c_x2, cord_y1, coord2);
			//coord1[0] = cord_x1; coord1[1] = cord_y1;
			//coord2[0] = c_x2+1; coord2[1] = cord_y1;
			minDist = dist2(coord1, coord2);
		}
		else {
			getLowerLeftPoint(cord_x1, cord_y1, coord1);
			getUpperRightPoint(c_x2, c_y2, coord2);
			//coord1[0] = cord_x1; coord1[1] = cord_y1;
			//coord2[0] = c_x2+1; coord2[1] = c_y2+1;
			minDist = dist2(coord1, coord2);				
		}

		getUpperRightPoint(cord_x2, cord_y2, coord1);
		getLowerLeftPoint(c_x1, c_y1, coord2);
		//coord1[0] = cord_x2+1; coord1[1] = cord_y2+1;
		//coord2[0] = c_x1; coord2[1] = c_y1;
		maxDist = dist2(coord1, coord2);

		return;
	}
}
//---------------------------------------------------------------------------
void Selection::DistanceF::getUpperRightPoint(unsigned x, unsigned y, double *coord) {

	coord[0] = x+1; //grid->fromUnitSquareToOriginalX(x+1);
	coord[1] = y+1; //grid->fromUnitSquareToOriginalY(y+1);
}
//---------------------------------------------------------------------------
void Selection::DistanceF::getUpperLeftPoint(unsigned x, unsigned y, double *coord) {

	coord[0] = x; 	 //grid->fromUnitSquareToOriginalX(x);
	coord[1] = y+1; //grid->fromUnitSquareToOriginalY(y+1);
}
//---------------------------------------------------------------------------
void Selection::DistanceF::getLowerLeftPoint(unsigned x, unsigned y, double *coord) {

	coord[0] = x;	//grid->fromUnitSquareToOriginalX(x);
	coord[1] = y;	//grid->fromUnitSquareToOriginalY(y);
}
//---------------------------------------------------------------------------
void Selection::DistanceF::getLowerRightPoint(unsigned x, unsigned y, double *coord) {

	coord[0] = x+1;	//grid->fromUnitSquareToOriginalX(x+1);
	coord[1] = y;	   //grid->fromUnitSquareToOriginalY(y);
}
//---------------------------------------------------------------------------
double Selection::DistanceF::dist2(double *c1, double *c2)
{
	double x = (c1[0] - c2[0]), y = (c1[1] - c2[1]);

	// return the square of the distance
	return (x*x*param1 + y*y*param2);
}
//---------------------------------------------------------------------------
void Selection::Result::ensureString(Selection *selection)
   // Ensure that a string is available
{
   //cout << "Ensuring string..." << endl;
	
	//cout << "flags is " << flags << endl;
	//cout << "?geo id is " << id << endl;	

   if (!(flags & stringAvailable)) {
		//cout << "Enter if" << endl;
      if (flags & idAvailable) {
         const char *start, *stop;
         if (~id) {
	    		if ( (selection->itter = selection->cache.find(id)) != selection->cache.end() ) {
					//cout << "First Case" << endl;
					value = selection->itter->second;
               flags |= typeAvailable;
					//cout << "value is " << value << endl;
					//cout << "flags is " << flags << endl; 
	    		}
	    		else if ( selection->runtime.getDatabase().getDictionary().lookupById(id, start, stop, type, subType) ) {
					//cout << "Second Case" << endl;
	         	value = string(start, stop);
		    		selection->cache.insert(pair<unsigned, string> (id, value));
	            flags |= typeAvailable;
					//cout << "value is " << value << endl;
					//cout << "flags is " << flags << endl;
	    		}
         } 
			else {
            value = "NULL";
         }
      } 
		else if (flags & booleanAvailable) {
         if (boolean)
            value = "true"; 
			else
            value = "false";
      } 
		else {
         value = "";
      }
      flags |= stringAvailable;
		//cout << "flags is " << flags << endl;
   }
}
//---------------------------------------------------------------------------
void Selection::Result::ensureType(Selection* selection)
   // Ensure that the type is available
{
   //std::cout << "Ensuring Type..." << std::endl;

   if (!(flags&typeAvailable)) {
      if (flags&idAvailable) {
         const char* start,*stop;
         if ((~id)&&(selection->runtime.getDatabase().getDictionary().lookupById(id,start,stop,type,subType))) {
            value=string(start,stop);
            flags|=stringAvailable;
         } else {
            type=Type::Literal; // XXX NULL type?
         }
      } else if (flags&booleanAvailable) {
         type=Type::Boolean;
      } else {
         type=Type::Literal;
      }
      flags|=typeAvailable;
   }
}
//---------------------------------------------------------------------------
void Selection::Result::ensureSubType(Selection* selection)
   // Ensure that the type is available
{

   //std::cout << "Ensuring SubType..." << std::endl;

   ensureType(selection);
   if (!(flags&subTypeAvailable)) {
      if ((type==Type::CustomLanguage)||(type==Type::CustomType)) {
         const char* start,*stop;
         Type::ID t; unsigned st;
         if (selection->runtime.getDatabase().getDictionary().lookupById(subType,start,stop,t,st)) {
            subTypeValue=string(start,stop);
         } else {
            subTypeValue.clear();
         }
      } else {
         subTypeValue.clear();
      }
      flags|=subTypeAvailable;
   }
}
//---------------------------------------------------------------------------
void Selection::Result::ensureBoolean(Selection *runtime)
   // Ensure that a boolean interpretation is available
{
   //cout << "Ensuring Boolean..." << endl;

   //cout << "flags is " << flags << endl;
	//cout << "flags & booleanAvailable is " << (flags & booleanAvailable) << endl;
	//cout << "!(flags & booleanAvailable) is " << !(flags & booleanAvailable) << endl;

   if (!(flags & booleanAvailable)) {
		//cout << "Enter If" << endl;
      ensureString(runtime);
      boolean = (value == "true");
      flags |= booleanAvailable;
   }
}
//---------------------------------------------------------------------------
void Selection::Result::setBoolean(bool v)
   // Set to a boolean value
{
   flags = booleanAvailable | typeAvailable;
   type = Type::Boolean;
   boolean = v;
}
//---------------------------------------------------------------------------
void Selection::Result::setId(unsigned v)
   // Set to an id value
{
   flags = idAvailable;
   id = v;
}
//---------------------------------------------------------------------------
void Selection::Result::setLiteral(const std::string& v)
   // Set to a string value
{
   flags=stringAvailable|typeAvailable;
   type=Type::Literal;
   value=v;
}
//---------------------------------------------------------------------------
void Selection::Result::setIRI(const std::string& v)
   // Set to a string value
{
   flags=stringAvailable|typeAvailable;
   type=Type::URI;
   value=v;
}
//---------------------------------------------------------------------------
Selection::Predicate::Predicate()
   : selection(0)
   // Constructor
{
}
//---------------------------------------------------------------------------
Selection::Predicate::~Predicate()
   // Destructor
{
}
//---------------------------------------------------------------------------
void Selection::Predicate::setSelection(Selection *s)
   // Set the selection
{
	//cout << "Predicate::setSelection" << endl;
   selection = s;
}
//---------------------------------------------------------------------------
bool Selection::Predicate::check()
   // Check the predicate
{
	//cout << "Check the predicate" << endl;
   Result r;
   eval(r);	 // evaluate the predicate
   r.ensureBoolean(selection);
   return r.boolean;
}
//---------------------------------------------------------------------------
Selection::BinaryPredicate::~BinaryPredicate()
   // Destructor
{
   delete left;
   delete right;
}
//---------------------------------------------------------------------------
void Selection::BinaryPredicate::setSelection(Selection* s)
   // Set the selection
{
   Predicate::setSelection(s);
   left->setSelection(s);
   right->setSelection(s);
}
//---------------------------------------------------------------------------
Selection::UnaryPredicate::~UnaryPredicate()
   // Destructor
{
   delete input;
}
//---------------------------------------------------------------------------
void Selection::UnaryPredicate::setSelection(Selection *s)
   // Set the selection
{
	//cout << "UnaryPredicate::setSelection" << endl;
   Predicate::setSelection(s);
   input->setSelection(s);
}
//---------------------------------------------------------------------------
void Selection::Or::eval(Result& result)
   // Evaluate the predicate
{
   result.setBoolean(left->check()||right->check());
}
//---------------------------------------------------------------------------
string Selection::Or::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "("+left->print(out)+")||("+right->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::And::eval(Result& result)
   // Evaluate the predicate
{
   result.setBoolean(left->check()&&right->check());
}
//---------------------------------------------------------------------------
string Selection::And::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "("+left->print(out)+")&&("+right->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::Equal::eval(Result& result)
   // Evaluate the predicate
{
   Result l,r;
   left->eval(l);
   right->eval(r);

   //add verification info
   selection->verified = selection->input->verified;

   //cout << "Equal selection verified: " << selection->verified << endl;

   // Cheap case first
   if (l.hasId()&&r.hasId()) {
      result.setBoolean(l.id==r.id);
	//if(l.id==r.id)
	//cout << "They are equal." << endl;
	//if(selection->verified){
	//cout << "They are verified." << endl;	
	//getchar();
	//}
      return;
   }

   //cout << "Selection::Equal." << endl;
   //getchar();

   // Now compare for real
   l.ensureString(selection);
   r.ensureString(selection);
   result.setBoolean(l.value==r.value);
}
//---------------------------------------------------------------------------
string Selection::Equal::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "("+left->print(out)+")==("+right->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::NotEqual::eval(Result& result)
   // Evaluate the predicate
{
   Result l,r;
   left->eval(l);
   right->eval(r);

   // Cheap case first
   if (l.hasId()&&r.hasId()) {
      result.setBoolean(l.id!=r.id);
      return;
   }

   // Now compare for real
   l.ensureString(selection);
   r.ensureString(selection);
   result.setBoolean(l.value!=r.value);
}
//---------------------------------------------------------------------------
string Selection::NotEqual::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "("+left->print(out)+")!=("+right->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::WithinF::getUpperRightPoint(unsigned x, unsigned y, double *coord) {
	coord[0] = (x+1)*(grid->cSide);
	coord[1] = (y+1)*(grid->cSide);
}
//---------------------------------------------------------------------------
void Selection::WithinF::getLowerLeftPoint(unsigned x, unsigned y, double *coord) {
	coord[0] = x*(grid->cSide);
	coord[1] = y*(grid->cSide);
}
//---------------------------------------------------------------------------
void Selection::WithinF::eval(Result &result)
	// Evaluate within predicate (filtering)
{
   Result i;
   input->eval(i);

   //std::cout << "Entered withinF selection..." << std::endl;
   //std::cout << "id: " << i.id << std::endl;

   // get c_id, gpos and level
   if (!getCellId(i.id, c_id, gpos, level)) {
	   //std::cout << "id: " << i.id << std::endl;
	   //cout << "WithinF: Not a spatial entity." << endl;
	   selection->verified = true;
	   result.setBoolean(false);	
	   return;
   }  
   
   if (!level) {  // it belongs to the whole grid
	   //std::cout << "id: " << i.id << std::endl;
	   //cout << "WithinF: Entity refers to the whole grid." << endl;
	   selection->verified = false;
	   result.setBoolean(true);	
	   return;
   }

   if (level < grid->MAX_LEVEL) {
		//std::cout << "id: " << i.id << std::endl;
		//cout << "WithinF: Not at max level." << endl;
		gpos -= (((unsigned) 1) << grid->b1);
		cX1 = grid->offsets[gpos][0];
		cY1 = grid->offsets[gpos][1];
		cX2 = grid->offsets[gpos][2];
		cY2 = grid->offsets[gpos][3];	
		//cout << "cX1: " << cX1 << " cY1: " << cY1 << " cX2: " << cX2 << " cY2: " << cY2 << endl;
   }
   else {
		//cout << "WithinF: At max level." << endl;
		cX1 = cX2 = grid->i2c[c_id].first;
		cY1 = cY2 = grid->i2c[c_id].second;
		//if(4104312831==i.id){
			//cout << "cX1: " << cX1 << " cY1: " << cY1 << " cX2: " << cX2 << " cY2: " << cY2 << endl;
			//getchar();
		//}
   }	

   // get lower left point coords in Unit Square
   getLowerLeftPoint(cX1, cY1, c1);   
   // get upper right point coords in Unit Square
   getUpperRightPoint(cX2, cY2, c2);

   // rectangle coordinates in the Unit Square   
   rx1 = ((Geo::Rectangle*) geo)->getP1()->getX(); 
   ry1 = ((Geo::Rectangle*) geo)->getP1()->getY(); 
   rx2 = ((Geo::Rectangle*) geo)->getP2()->getX(); 
   ry2 = ((Geo::Rectangle*) geo)->getP2()->getY();

   if (c1[0] >= rx1 && c2[0] <= rx2 && c1[1] >= ry1 && c2[1] <= ry2) {  // cell falls inside the given window
		selection->verified = true;
		result.setBoolean(true);
   	//cout << "rx1: " << rx1 << " ry1: " << ry1 << " rx2: " << rx2 << " ry2: " << ry2 << endl;
		//getchar();
		//cout << "ID: " << i.id << endl;

		//if(4104314751==i.id){
   		//cout << "c1[0]: " << c1[0] << " c1[1]: " << c1[1] << " c2[0]: " << c2[0] << " c2[1]: " << c2[1] << endl;
			//cout << "rx1: " << rx1 << " ry1: " << ry1 << " rx2: " << rx2 << " ry2: " << ry2 << endl;
			//getchar();
		//}
		return;
   }
   
   if (c1[0] > rx2 || c2[0] < rx1 || c1[1] > ry2 || c2[1] < ry1) {  // cell is out of the window
	   selection->verified = true;
	   //cout << "WithinF2 res: " << 0 << " ver: " << selection->verified << endl;
	   selection->filtered++;
	   result.setBoolean(false);
	   //getchar();
	   return;
   }

   // cell overlaps with the given window
   selection->verified = false;
   result.setBoolean(true);
   //cout << "WithinF3 res: " << res << " ver: " << selection->verified << endl;
}
//---------------------------------------------------------------------------
void Selection::DistanceF::eval(Result &result)
	// Evaluate distance predicate (filtering)
{
	Result r1,r2;
	arg1->eval(r1);
	arg2->eval(r2);

	/*
        bool sf = false;
	if(r1.id == r2.id){
		cout << "They are the same." << endl;
		sf = true;
	}
	*/

   	//std::cout << "Entered DistanceF selection..." << std::endl;

	// get the cell id, g_pos and the level of the grid - if at least one of the entities is not spatial return immediately
	if (!getCellId(r1.id, lcID, g_pos1, level1)) { result.setBoolean(false); return; }
	if (!getCellId(r2.id, rcID, g_pos2, level2)) { result.setBoolean(false); return; }

	if (!level1 || !level2) {  // at least one entity refers to the whole grid
		selection->verified = false;
		result.setBoolean(true);
		return;
	}

	// if both cells refer to MAX_LEVEL
	if ((level1 == grid->MAX_LEVEL) && (level2 == grid->MAX_LEVEL)) {

		//cout << "Both cells refer to MAX_LEVEL." << endl;

		// get coordinates (lower left points)
		leftCellX1 = leftCellX2 = grid->i2c[lcID].first;
		leftCellY1 = leftCellY2 = grid->i2c[lcID].second;
		rightCellX1 = rightCellX2 = grid->i2c[rcID].first; 
		rightCellY1 = rightCellY2 = grid->i2c[rcID].second;

		getMinMaxDist(leftCellX1, leftCellY1, leftCellX2, leftCellY2, rightCellX1, rightCellY1, rightCellX2, rightCellY2, minDist, maxDist);
		res = evalPredicate(minDist, maxDist, selection->verified);
		result.setBoolean(res);
		if (!res) selection->filtered++; 
		//cout << "Edw1." << endl;

		/*
		if(sf){

			cout << "REs0: " << res << endl;
			getchar();
		}
		*/

		return;
	}
	else if (level1 == grid->MAX_LEVEL) {  // left cell refers to the max level

		//cout << "Left cell refers to MAX_LEVEL." << endl;

		// lower left offsets for the left cell
		leftCellX1 = leftCellX2 = grid->i2c[lcID].first;
		leftCellY1 = leftCellY2 = grid->i2c[lcID].second;

		g_pos2 -= (((unsigned) 1) << grid->b1);
		// lower left and upper right offsets for the right cell
		rightCellX1 = grid->offsets[g_pos2][0]; 
		rightCellY1 = grid->offsets[g_pos2][1]; 
		rightCellX2 = grid->offsets[g_pos2][2]; 
		rightCellY2 = grid->offsets[g_pos2][3];

		getMinMaxDist(leftCellX1, leftCellY1, leftCellX2, leftCellY2, rightCellX1, rightCellY1, rightCellX2, rightCellY2, minDist, maxDist);
		res = evalPredicate(minDist, maxDist, selection->verified);
		result.setBoolean(res);
		if (!res) selection->filtered++; 
		//cout << "Edw2." << endl;

		/*
		if(sf){

			cout << "REs0: " << res << endl;
			getchar();
		}
		*/

		return;
	}
	else if (level2 == grid->MAX_LEVEL) {  // right cell refers to the max level

		//cout << "Right cell refers to MAX_LEVEL." << endl;

		// lower left offsets for the right cell
		rightCellX1 = rightCellX2 = grid->i2c[rcID].first; 
		rightCellY1 = rightCellY2 = grid->i2c[rcID].second;

		g_pos1 -= (((unsigned) 1) << grid->b1);
		// lower left and upper right offsets for the left cell
		leftCellX1 = grid->offsets[g_pos1][0], 
		leftCellY1 = grid->offsets[g_pos1][1], 
		leftCellX2 = grid->offsets[g_pos1][2], 
		leftCellY2 = grid->offsets[g_pos1][3];

		getMinMaxDist(leftCellX1, leftCellY1, leftCellX2, leftCellY2, rightCellX1, rightCellY1, rightCellX2, rightCellY2, minDist, maxDist);
		res = evalPredicate(minDist, maxDist, selection->verified);
		result.setBoolean(res);
		if (!res) selection->filtered++; 		
		//cout << "Edw3." << endl;

		/*
		if(sf){

			cout << "REs0: " << res << endl;
			getchar();
		}
		*/

		return;
	}
	else {  // both cells refer to a higher level

		//cout << "Both cells refer to a higher level." << endl;

		g_pos1 -= (((unsigned) 1) << grid->b1);
		// lower left and upper right offsets for the left cell
		leftCellX1 = grid->offsets[g_pos1][0];
		leftCellY1 = grid->offsets[g_pos1][1]; 
		leftCellX2 = grid->offsets[g_pos1][2];
		leftCellY2 = grid->offsets[g_pos1][3];

		g_pos2 -= (((unsigned) 1) << grid->b1);
		// lower left and upper right offsets for the right cell
		rightCellX1 = grid->offsets[g_pos2][0]; 
		rightCellY1 = grid->offsets[g_pos2][1]; 
		rightCellX2 = grid->offsets[g_pos2][2];
		rightCellY2 = grid->offsets[g_pos2][3];

		getMinMaxDist(leftCellX1, leftCellY1, leftCellX2, leftCellY2, rightCellX1, rightCellY1, rightCellX2, rightCellY2, minDist, maxDist);
		res = evalPredicate(minDist, maxDist, selection->verified);
		result.setBoolean(res);
		if (!res) selection->filtered++; 
		//cout << "Edw4." << endl;

		/*
		if(sf){

			cout << "REs0: " << res << endl;
			getchar();
		}
		*/

		return;
	}
}
//---------------------------------------------------------------------------
Selection::DistanceF::~DistanceF()
   // Destructor
{
   delete arg1;
   delete arg2;
}
//---------------------------------------------------------------------------
void Selection::DistanceF::setSelection(Selection *s)
   // Set the selection
{
   Predicate::setSelection(s);
   arg1->setSelection(s);
   arg2->setSelection(s);
   arg3->setSelection(s);
}
#ifdef HASGEO
void Selection::Within::eval(Result &result)
	// Evaluate within predicate
{
	//cout << "Evaluate WITHIN predicate" << endl;
#if GEOSTORE
	cout << "Den mpainei edw" << endl;
   //get the id
   Result r;
   input->eval(r);

   cnt++;

   if((r.id < gr->CBOT) && gr->grid[r.id]){//it is a filter on pos id

		//get the normalized coordinates
		minX = gr->i2c[r.id].first*gr->cSide;
		maxX = minX + gr->cSide;
		minY = gr->i2c[r.id].second*gr->cSide;
		maxY = minY + gr->cSide;

		/*
		if(r.id==10674849){
			cout << "X1: " << minX  << " Y1: " << minY << " X2: " << maxX << " Y2: " << maxY << endl;
			cout << "X1: " << ((Geo::Rectangle*)geo2)->getP1()->getX()  << " Y1: " << ((Geo::Rectangle*)geo2)->getP1()->getY() << " X2: " << ((Geo::Rectangle*)geo2)->getP2()->getX() << " Y2: " << ((Geo::Rectangle*)geo2)->getP2()->getY() << endl;
		}
		*/
		//getchar();

	   //check WITHIN predicate
	   if(minX>=((Geo::Rectangle*)geo2)->getP1()->getX() && maxX<=((Geo::Rectangle*)geo2)->getP2()->getX() && minY>=((Geo::Rectangle*)geo2)->getP1()->getY() && maxY<=((Geo::Rectangle*)geo2)->getP2()->getY()){

			result.setBoolean(true);
			//std::cout << "res: " << 1 << std::endl;
			//exit(-1);

			//std::cout << "long: " << n1 << " lat: " << n2 << std::endl;

			//deallocate
			//delete[] points;

			selection->verified = true;
			return;
	   }
	   else if(maxX<((Geo::Rectangle*)geo2)->getP1()->getX() || minX>((Geo::Rectangle*)geo2)->getP2()->getX() || maxY<((Geo::Rectangle*)geo2)->getP1()->getY() || minY>((Geo::Rectangle*)geo2)->getP2()->getY()){

			result.setBoolean(false);
			selection->verified = true;
			return;
	   }	 
		
	   //if(selection->input->verified){

		//cout << "Selection: This should not happen." << endl;
		//std::cout << "long: " << n1 << " lat: " << n2 << std::endl;
	   //}

	   selection->verified = false;
	   //std::cout << "res: " << 0 << std::endl;
	   result.setBoolean(true); 
	   return;	
   }
   else{

	   //if(r.id==10674849){
	   //		cout << "Here!" << endl;
	   //}

	   //if the output of the previous operator is verified
	   if(selection->input->verified){

		//std::cout << "Within result is verified!" << std::endl;

		result.setBoolean(true);
		selection->verif++;

		return;		
	   }

	   //std::cout << "Getting coordinate." << std::endl;	
	   //retrieve the geometry
	     r.ensureString(selection);

	   //get the coordinates (original coordinates + 180)
	   getLongLat(r.value,points,size);

	   //get the MBR of the geometry
	   getRegion(points,size,minX,maxX,minY,maxY);

	   //std::cout << "Entered within selection..." << std::endl;
	   //std::cout << "long: " << n1 << " lat: " << n2 << std::endl;
	   //getchar();

	   //check WITHIN predicate
	   if(minX>=((Geo::Rectangle*)geo)->getP1()->getX() && maxX<=((Geo::Rectangle*)geo)->getP2()->getX() && minY>=((Geo::Rectangle*)geo)->getP1()->getY() && maxY<=((Geo::Rectangle*)geo)->getP2()->getY()){

			result.setBoolean(true);
			//std::cout << "res: " << 1 << std::endl;
			//exit(-1);

			//std::cout << "long: " << n1 << " lat: " << n2 << std::endl;

			//deallocate
			//delete[] points;

			return;
	   }
	 
	   //if(selection->input->verified){

		//cout << "Selection: This should not happen." << endl;
		//std::cout << "long: " << n1 << " lat: " << n2 << std::endl;
	   //}

	   //std::cout << "res: " << 0 << std::endl;
	   result.setBoolean(false); 
   }
#else

   // get the id
   Result r;
   input->eval(r);

   cnt++;

   // if the output of the previous operator is verified
   if (selection->input->verified) {
		//cout << "Pass Within WITHOUT geometry retrieval" << endl;
		result.setBoolean(true);
		selection->verif++;
		return;		
   }

   //std::cout << "Getting coordinates" << std::endl;	

   // retrieve the geometry (get the ?geo string::value)
   r.ensureString(selection);

   // get the coordinates (original coordinates + 180)
   getLongLat(r.value, points, size);

   // get the MBR of the geometry
   getRegion(points, size, minX, maxX, minY, maxY);

   //std::cout << "Entered within selection..." << std::endl;
   //std::cout << "long: " << n1 << " lat: " << n2 << std::endl;
   //getchar();

   // check WITHIN predicate
   if (minX >= ((Geo::Rectangle*) geo)->getP1()->getX() && maxX <= ((Geo::Rectangle *)geo)->getP2()->getX() && 
       minY >= ((Geo::Rectangle*) geo)->getP1()->getY() && maxY <= ((Geo::Rectangle *)geo)->getP2()->getY()) {
 
		result.setBoolean(true);
      //cout << "Pass Within WITH geometry retrieval" << endl;

		//std::cout << "res: " << 1 << std::endl;
		//exit(-1);

		//std::cout << "long: " << n1 << " lat: " << n2 << std::endl;

		//deallocate
		//delete[] points;

		return;
   }
 
   //if(selection->input->verified){

	//cout << "Selection: This should not happen." << endl;
	//std::cout << "long: " << n1 << " lat: " << n2 << std::endl;
   //}

   //std::cout << "res: " << 0 << std::endl;

   result.setBoolean(false); 
   //cout << "Do not pass Within" << endl;   

   //deallocate
   //delete[] points;
#endif
}
//---------------------------------------------------------------------------
void Selection::Within::getRegion(double *points, unsigned size, double &minX, double &maxX, double &minY, double &maxY)
   // Get the coordinates of the MBR for a given geometry
{
	minX = minY = 1000000;	//just a big value
	maxX = maxY = 0;	//the smallest possible, given that the coordinates are normalized in [0,360]

	//cout << "Size: " << size << endl;

	for(unsigned i=0;i<size;i++){

		if(i%2){//latitude

			if(minY>points[i])
				minY = points[i];

			if(maxY<points[i])
				maxY = points[i];
		}
		else{//longitude

			if(minX>points[i])
				minX = points[i];

			if(maxX<points[i])
				maxX = points[i];	
		}
	}

	//cout << "minX: " << minX << " maxX: " << maxX << " minY: " << minY << " maxY: " << maxY << endl;
	//getchar();
}
//---------------------------------------------------------------------------
string Selection::Within::print(PlanPrinter &out)
   // Print the predicate (debugging only)
{
   return "Within("+input->print(out)+") IN "+geo->toString();
}
//---------------------------------------------------------------------------
unsigned Selection::Within::countPoints(const string &geom, const string &dlm)
{
	unsigned num = 0;

	size_t pos = geom.find(dlm);

	while(pos!=std::string::npos){

		pos = geom.find(dlm,pos+1);
		num++;
	}
	
	if(num)	num++;
	
	//cout << "Found " << num << " points in geometry: " << geom << endl;

	return num;
}
//---------------------------------------------------------------------------
void Selection::Within::getLongLat(const std::string &geometry, double* &points, unsigned &size)
{
	size = (geometry.length() - sizeof(unsigned))/sizeof(double);

	const char *data = geometry.data();
	
	//points = new double[size];
	
	points = (double*) &data[sizeof(unsigned)];

	//ignore the geometry type
	//memcpy(points,&data[sizeof(unsigned)],size*sizeof(double));
	
/*
	string cm = ",", dlm = "|", endS = ")", 
	       point = "POINT (", 
	       line = "LINESTRING (", 
	       polygon = "POLYGON (", 
	       multipoint = "MULTIPOINT (";

	unsigned pLength = point.length(), lLength = line.length(), plLength = polygon.length(), mpLength = multipoint.length();

	//cout << "Input Geometry: " << str << endl;

	bool p=false,l=false,pl=false,mp=false;

	//size_t pos = geometry.find(point);
	
	//if (pos==std::string::npos){
	if(geometry[3]!='N'){

		//cout << "It is not a point." << endl;

		//pos = geometry.find(line);

		//if (pos==std::string::npos){
		if(geometry[3]!='E'){

			//cout << "It is not a line." << endl;

			//pos = geometry.find(polygon);

			//if (pos==std::string::npos){
			if(geometry[3]!='Y'){

				//cout << "It is not a polygon." << endl;

				//pos = geometry.find(multipoint);

				//if (pos==std::string::npos){	
				if(geometry[3]!='T'){

					cout << "Could not recognize geometry." << endl;
					cout << geometry << endl;
					exit(-1);
				}

				mp = true;
			}
			else	pl = true;
		}
		else l = true;
	}
	else p = true;
	
	string lon,lat;

	if(p){//point, only two coords

		size = 2;
		points = new double[size];

		size_t pos = geometry.find(cm,pLength);

		//cout << "pos: " << pos << endl;

		lon = geometry.substr(pLength,pos-pLength);
		lat = geometry.substr (pos+1,geometry.find(endS,pos+1)-pos-1);

		//cout << "POINT: " << geometry << " has Longitude: " << lon << " and Latitude: " << lat << endl;

		stringstream ss1(stringstream::in|stringstream::out),ss2(stringstream::in|stringstream::out);

		ss1 << lon;
		ss1 >> points[0];
		ss2 << lat;
		ss2 >> points[1];

		//cout << "POINT: " << geometry << " has Longitude: " << points[0] << " and Latitude: " << points[1] << endl;
		//getchar();
	}
	else if(l){//line string, many points

		size = 2*countPoints(geometry,dlm); 
		points = new double[size];
	
		//first comma
		size_t pos = geometry.find(cm,lLength);
		size_t pos2 = geometry.find(dlm,lLength);

		//get first point
		lon = geometry.substr(lLength,pos-lLength);
		lat = geometry.substr(pos+1,pos2-pos-1);
		
		stringstream s1(stringstream::in|stringstream::out),s2(stringstream::in|stringstream::out);
		s1 << lon;
		s1 >> points[0];
		s2 << lat;
		s2 >> points[1];

		//cout << "First point of LINE " << geometry << " has Longitude: " << points[0] << " and Latitute: " << points[1] << endl; 

		//first delimeter
		pos = geometry.find(dlm,lLength);
		
		//second comma
		pos2 = geometry.find(cm,pos+1);
		unsigned ind = 2;
		size_t temp;	
	
		//get remaining points
		while(geometry.find(dlm,pos2+1)!=std::string::npos){

			lon = geometry.substr(pos+1,pos2-pos-1);
			pos = geometry.find(dlm,pos2+1);
			lat = geometry.substr (pos2+1,pos-pos2-1);

			stringstream ss1(stringstream::in|stringstream::out),ss2(stringstream::in|stringstream::out);
			ss1 << lon;
			ss1 >> points[ind];
			ind++;
			ss2 << lat;
			ss2 >> points[ind];
			ind++;		
			
			//next comma
			pos2 = geometry.find(cm,pos+1);	
		}

		//get last point
		pos2 = geometry.find(cm,pos+1);
		
		lon = geometry.substr(pos+1,pos2-pos+1);
		lat = geometry.substr (pos2+1,geometry.find(endS,pos2+1)-pos2-1);

		stringstream sss1(stringstream::in|stringstream::out),sss2(stringstream::in|stringstream::out);
		sss1 << lon;
		sss1 >> points[ind];
		ind++;
		sss2 << lat;
		sss2 >> points[ind];
		ind++;

		//cout << "Last point of LINE " << geometry << " has Longitude: " << points[ind-2] << " and Latitute: " << points[ind-1] << endl;

		if(ind!=size){

			cout << "Something is wrong with the line index." << endl;
			cout << "Index: " << ind << " Size: " << size << endl;
			exit(-1);
		}	

		//getchar();

	}
	else if(pl){//polygon, many points
	
		size = 2*countPoints(geometry,dlm); 
		points = new double[size];
	
		//first comma
		size_t pos = geometry.find(cm,plLength);
		size_t pos2 = geometry.find(dlm,plLength);

		//get first point
		lon = geometry.substr(plLength,pos-plLength);
		lat = geometry.substr(pos+1,pos2-pos-1);
		
		stringstream s1(stringstream::in|stringstream::out),s2(stringstream::in|stringstream::out);
		s1 << lon;
		s1 >> points[0];
		s2 << lat;
		s2 >> points[1];

		//cout << "First point of POLYGON " << geometry << " has Longitude: " << points[0] << " and Latitute: " << points[1] << endl;

		//first delimeter
		pos = geometry.find(dlm,plLength);
		
		//second comma
		pos2 = geometry.find(cm,pos+1); 
		unsigned ind = 2;
		size_t temp;	

		//get remaining points
		while(geometry.find(dlm,pos2+1)!=std::string::npos){

			lon = geometry.substr(pos+1,pos2-pos-1);

			//next delimeter
			pos = geometry.find(dlm,pos2+1);

			lat = geometry.substr (pos2+1,pos-pos2-1);

			//cout << "Lon: " << lon << " Lat: " << lat << endl;

			stringstream ss1(stringstream::in|stringstream::out),ss2(stringstream::in|stringstream::out);
			ss1 << lon;
			ss1 >> points[ind];
			ind++;
			ss2 << lat;
			ss2 >> points[ind];
			ind++;		
			
			//next comma
			pos2 = geometry.find(cm,pos+1);	
		}

		//cout << "Index: " << ind << endl;

		//get last point
		pos2 = geometry.find(cm,pos+1);
		
		lon = geometry.substr(pos+1,pos2-pos-1);
		lat = geometry.substr (pos2+1,geometry.find(endS,pos2+1)-pos2-1);

		stringstream sss1(stringstream::in|stringstream::out),sss2(stringstream::in|stringstream::out);
		sss1 << lon;
		sss1 >> points[ind];
		sss2 << lat;
		ind++;
		sss2 >> points[ind];
		ind++;

		//cout << "Last point of POLYGON " << geometry << " has Longitude: " << points[ind-2] << " and Latitute: " << points[ind-1] << endl;

		if(ind!=size){

			cout << "Something is wrong with the polygon index." << endl;
			cout << "Index: " << ind << " Size: " << size << endl;
			exit(-1);
		}	
	}
	else{//multipoint, many points
	
		size = 2*countPoints(geometry,dlm); 
		points = new double[size];
	
		//first comma
		size_t pos = geometry.find(cm,mpLength);
		size_t pos2 = geometry.find(dlm,mpLength);

		//get first point
		lon = geometry.substr(mpLength,pos-mpLength);
		lat = geometry.substr(pos+1,pos2-pos-1);
		
		stringstream s1(stringstream::in|stringstream::out),s2(stringstream::in|stringstream::out);
		s1 << lon;
		s1 >> points[0];
		s2 << lat;
		s2 >> points[1];

		//cout << "First point of MULTIPOINT " << geometry << " has Longitude: " << points[0] << " and Latitute: " << points[1] << endl;

		//first delimeter
		pos = geometry.find(dlm,mpLength);
		
		//second comma
		pos2 = geometry.find(cm,pos+1); 
		unsigned ind = 2;
		size_t temp;	

		//get remaining points
		while(geometry.find(dlm,pos2+1)!=std::string::npos){

			lon = geometry.substr(pos+1,pos2-pos-1);

			//next delimeter
			pos = geometry.find(dlm,pos2+1);

			lat = geometry.substr (pos2+1,pos-pos2-1);

			//cout << "Lon: " << lon << " Lat: " << lat << endl;

			stringstream ss1(stringstream::in|stringstream::out),ss2(stringstream::in|stringstream::out);
			ss1 << lon;
			ss1 >> points[ind];
			ind++;
			ss2 << lat;
			ss2 >> points[ind];
			ind++;		
			
			//next comma
			pos2 = geometry.find(cm,pos+1);	
		}

		//cout << "Index: " << ind << endl;

		//get last point
		pos2 = geometry.find(cm,pos+1);
		
		lon = geometry.substr(pos+1,pos2-pos-1);
		lat = geometry.substr (pos2+1,geometry.find(endS,pos2+1)-pos2-1);

		stringstream sss1(stringstream::in|stringstream::out),sss2(stringstream::in|stringstream::out);
		sss1 << lon;
		sss1 >> points[ind];
		sss2 << lat;
		ind++;
		sss2 >> points[ind];
		ind++;

		//cout << "Last point of MULTIPOINT " << geometry << " has Longitude: " << points[ind-2] << " and Latitute: " << points[ind-1] << endl;

		if(ind!=size){

			cout << "Something is wrong with the multipoint index." << endl;
			cout << "Index: " << ind << " Size: " << size << endl;
			exit(-1);
		}	
	}
*/

}
//---------------------------------------------------------------------------
double Selection::Distance::dot(double* A, double*B)
	//Compute the dot product of vectors A and B
{
        return A[0] * B[0] + A[1] * B[1];
}
//---------------------------------------------------------------------------
double Selection::Distance::distance2(double* A, double* B)
	//Compute the square of the distance from A to B
{
        double d1 = A[0] - B[0];
        double d2 = A[1] - B[1];
        return d1*d1+d2*d2;
}
//---------------------------------------------------------------------------
double Selection::Distance::pointLineDist2(double* point, double* segment)
    	//Compute the square of the minimum distance between a point and a line segment
{       
	double v[2],w[2];

	//Vector v = S.P1 - S.P0;
	v[0] = segment[2] - segment[0];
	v[1] = segment[3] - segment[1];

	//Vector w = P - S.P0;
	w[0] = point[0] - segment[0];
	w[1] = point[1] - segment[1];

	double c1 = dot(w,v);
     	if ( c1 <= 0 )
          	return distance2(point,segment);	//distance2(P, S.P0);

     	double c2 = dot(v,v);
     	if ( c2 <= c1 )
          	return distance2(point,&segment[2]);	//distance2(P, S.P1);

     	double b = c1 / c2;
	double pb[2];

	//Point Pb = S.P0 + b * v;
	pb[0] = segment[0] + b*v[0];
	pb[1] = segment[1] + b*v[1];

    	return distance2(point, pb);	//return distance2(P, Pb);
}
//---------------------------------------------------------------------------
double Selection::Distance::pointLineStringDist2(double* point, double* line, unsigned size)
    	//Compute the square of the minimum distance between a point and a line segment
{       
	double d, mD = 100000000;	//the minimum distance
	
	size -= 2;

	//for each segment of the given linestring
	for(unsigned i=0;i<size;i+=2){

		//distance between the line and the point
		d = pointLineDist2(point,&line[i]);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::Distance::lineStringLineStringDist2(double* line1, unsigned size1, double* line2, unsigned size2)
    	//Compute the square of the minimum distance between a point and a line segment
{       
	double d, mD = 100000000;	//the minimum distance
	
	size1 -= 2;
	size2 -= 2;

	//for each segment of the first linestring
	for(unsigned i=0;i<size1;i+=2){

		//for each segment of the second linestring
		for(unsigned j=0;j<size2;j+=2){

			//distance between the lines
			d = lineLineDist2(&line1[i],&line2[j]);

			if(!d) return 0;

			if(d<mD) mD = d;
		}
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::Distance::lineStringPolygonDist2(double* line, unsigned size1, double* polygon, unsigned size2)
    	//Compute the square of the minimum distance between a point and a line segment
{       
	double d, mD = 100000000;	//the minimum distance
	
	size1 -= 2;

	//for each segment of the given linestring
	for(unsigned i=0;i<size1;i+=2){

		d = linePolygonDist2(&line[i],polygon,size2);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::Distance::lineLineDist2(double* seg1, double* seg2)
    	//Compute the square of the minimum distance between two line segments
{ 
	double u[2],v[2],w[2];

	//Vector   u = S1.P1 - S1.P0;
	u[0] = seg1[2] - seg1[0];
	u[1] = seg1[3] - seg1[1];

	//Vector   v = S2.P1 - S2.P0;
	v[0] = seg2[2] - seg2[0];
	v[1] = seg2[3] - seg2[1];
	
	//Vector   w = S1.P0 - S2.P0;
	w[0] = seg1[0] - seg2[0];
	w[1] = seg1[1] - seg2[1];

	double    a = dot(u,u);         // always >= 0
	double    b = dot(u,v);
	double    c = dot(v,v);         // always >= 0
	double    d = dot(u,w);
	double    e = dot(v,w);
	double    D = a*c - b*b;        // always >= 0
	double    sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
	double    tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

	// compute the line parameters of the two closest points
	if (D < 0.00000001) { 		// the lines are almost parallel
		sN = 0.0;         	// force using point P0 on segment S1
		sD = 1.0;         	// to prevent possible division by 0.0 later
		tN = e;
		tD = c;
	}
	else {                 		// get the closest points on the infinite lines
		sN = (b*e - c*d);
		tN = (a*e - b*d);
		
		if (sN < 0.0) {        	// sc < 0 => the s=0 edge is visible
	    		sN = 0.0;
	    		tN = e;
	    		tD = c;
		}
		else if (sN > sD) {	// sc > 1  => the s=1 edge is visible
	    		sN = sD;
	    		tN = e + b;
	    		tD = c;
		}
	}

	if (tN < 0.0) {           	// tc < 0 => the t=0 edge is visible

		tN = 0.0;

		// recompute sc for this edge
		if (-d < 0.0)
	    		sN = 0.0;
		else if (-d > a)
	    		sN = sD;
		else {
	    		sN = -d;
	    		sD = a;
		}
	}
	else if (tN > tD) {      	// tc > 1  => the t=1 edge is visible

		tN = tD;
	
		// recompute sc for this edge
		if ((-d + b) < 0.0)
	    		sN = 0;
		else if ((-d + b) > a)
		    sN = sD;
		else {
	    		sN = (-d +  b);
	    		sD = a;
		}
	}

	// finally do the division to get sc and tc
	sc = (abs(sN) < 0.00000001 ? 0.0 : sN / sD);
	tc = (abs(tN) < 0.00000001 ? 0.0 : tN / tD);

	// get the difference of the two closest points
	double dp[2];

	//Vector   dP = w + (sc * u) - (tc * v);  // =  S1(sc) - S2(tc)
	dp[0] = w[0] + sc *u[0] - tc*v[0];
	dp[1] = w[1] + sc *u[1] - tc*v[1];

	return dot(dp,dp);   // return the square of the minimum distance
}
//---------------------------------------------------------------------------
double Selection::Distance::linePolygonDist2(double* line, double* polygon, unsigned size)
    	//Compute the square of the minimum distance between a line segment and a polygon
{ 
	double d, mD = 100000000;	//the minimum distance
	
	size -= 2;

	//for each line of the given polygon
	for(unsigned i=0;i<size;i+=2){

		//distance between the lines
		d = lineLineDist2(line,&polygon[i]);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::Distance::pointPolygonDist2(double* point, double* polygon, unsigned size)
    	//Compute the square of the minimum distance between a line segment and a polygon
{ 
	double d, mD = 100000000;	//the minimum distance
	
	size -= 2;

	//for each line of the given polygon
	for(unsigned i=0;i<size;i+=2){

		//distance between the point and the line
		d = pointLineDist2(point,&polygon[i]);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::Distance::polygonPolygonDist2(double* polygon1, double* polygon2, unsigned size1, unsigned size2)
    	//Compute the square of the minimum distance between a line segment and a polygon
{ 
	double d, mD = 100000000;	//the minimum distance
	
	size1 -= 2;
	size2 -= 2;

	//for each line of the first polygon
	for(unsigned i=0;i<size1;i+=2){

		//for each line of the second polygon
		for(unsigned j=0;j<size2;j+=2){

			//distance between the lines
			d = lineLineDist2(&polygon1[i],&polygon2[j]);

			if(!d) return 0;

			if(d<mD) mD = d;
		}
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::Distance::pointMultipointDist2(double* point, double* multipoint, unsigned size)
    	//Compute the square of the minimum distance between a line segment and a polygon
{ 
	double d, mD = 100000000;	//the minimum distance
	
	//for each point of the given multipoint
	for(unsigned i=0;i<size;i+=2){

		//distance between the points
		d = distance2(point,&multipoint[i]);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::Distance::lineStringMultipointDist2(double* line, unsigned size1, double* multipoint, unsigned size2)
    	//Compute the square of the minimum distance between a line segment and a polygon
{ 
	double d, mD = 100000000;	//the minimum distance
	
	size1 -= 2;

	//for each segment of the given linestring
	for(unsigned i=0;i<size1;i+=2){

		//distance between the point and the the line
		d = lineMultipointDist2(&line[i],multipoint,size2);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::Distance::lineMultipointDist2(double* line, double* multipoint, unsigned size)
    	//Compute the square of the minimum distance between a line segment and a polygon
{ 
	double d, mD = 100000000;	//the minimum distance
	
	//for each point of the given multipoint
	for(unsigned i=0;i<size;i+=2){

		//distance between the point and the the line
		d = pointLineDist2(&multipoint[i],line);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::Distance::polygonMultipointDist2(double* polygon, double* multipoint, unsigned size1, unsigned size2)
    	//Compute the square of the minimum distance between a line segment and a polygon
{ 
	double d, mD = 100000000;	//the minimum distance
	
	size1 -= 2;

	//for each line of the given polygon
	for(unsigned i=0;i<size1;i+=2){

		//for each point of the given multipoint
		for(unsigned j=0;j<size2;j+=2){

			//distance between the point and the the line
			d = pointLineDist2(&multipoint[j],&polygon[i]);

			if(!d) return 0;

			if(d<mD) mD = d;
		}
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::Distance::multipointMultipointDist2(double* multipoint1, double* multipoint2, unsigned size1, unsigned size2)
    	//Compute the square of the minimum distance between a line segment and a polygon
{ 
	double d, mD = 100000000;	//the minimum distance
	
	//for each point of the first multipoint
	for(unsigned i=0;i<size1;i+=2){

		//for each point of the second multipoint
		for(unsigned j=0;j<size2;j+=2){

			//distance between the points
			d = distance2(&multipoint1[i],&multipoint2[j]);

			if(!d) return 0;

			if(d<mD) mD = d;
		}
	}

	return mD;
}
//---------------------------------------------------------------------------
unsigned Selection::Distance::countPoints(const string& geom, const string& dlm)
{
	unsigned num = 0;

	size_t pos = geom.find(dlm);

	while(pos!=std::string::npos){

		pos = geom.find(dlm,pos+1);
		num++;
	}
	
	if(num)	num++;
	
	//cout << "Found " << num << " points in geometry: " << geom << endl;

	return num;
}
//---------------------------------------------------------------------------
unsigned Selection::Distance::getLongLat(const std::string &geometry, double* &points, unsigned &size)
	// Return the type of geometry (0: point, 1:line, 2:polygon, 3:multipoint)
{
	size = (geometry.length() - sizeof(unsigned))/sizeof(double);

	if(size%2){
		cout << "Selection: Sth is wrong with the size of the geometry. Size is " << size << endl;
		exit(-1);
	}

	const char *data = geometry.data();
	
	//points = new double[size];

	points = (double*) &data[sizeof(unsigned)];

	//ignore the geometry type
	//memcpy(points,&data[sizeof(unsigned)],size*sizeof(double));

	return *(unsigned*) data; 
	
/*
	string cm = ",", dlm = "|", endS = ")", 
	       point = "POINT (", 
	       line = "LINESTRING (", 
	       polygon = "POLYGON (", 
	       multipoint = "MULTIPOINT (";

	unsigned pLength = point.length(), lLength = line.length(), plLength = polygon.length(), mpLength = multipoint.length();

	//cout << "Input Geometry: " << str << endl;

	bool p=false,l=false,pl=false,mp=false;

	size_t pos = geometry.find(point);
	
	if (pos==std::string::npos){

		//cout << "It is not a point." << endl;

		pos = geometry.find(line);

		if (pos==std::string::npos){

			//cout << "It is not a line." << endl;

			pos = geometry.find(polygon);

			if (pos==std::string::npos){

				//cout << "It is not a polygon." << endl;

				pos = geometry.find(multipoint);

				if (pos==std::string::npos){

					cout << "Could not recognize geometry." << endl;
					cout << geometry << endl;
					exit(-1);
				}

				mp = true;
			}
			else	pl = true;
		}
		else l = true;
	}
	else p = true;
	
	string lon,lat;

	if(p){//point, only two coords

		size = 2;
		points = new double[size];

		pos = geometry.find(cm,pLength);

		//cout << "pos: " << pos << endl;

		lon = geometry.substr(pLength,pos-pLength);
		lat = geometry.substr (pos+1,geometry.find(endS,pos+1)-pos-1);

		//cout << "POINT: " << geometry << " has Longitude: " << lon << " and Latitude: " << lat << endl;

		stringstream ss1(stringstream::in|stringstream::out),ss2(stringstream::in|stringstream::out);

		ss1 << lon;
		ss1 >> points[0];
		ss2 << lat;
		ss2 >> points[1];

		return 0;

		//cout << "POINT: " << geometry << " has Longitude: " << points[0] << " and Latitude: " << points[1] << endl;
		//getchar();
	}
	else if(l){//line string, many points

		size = 2*countPoints(geometry,dlm); 
		points = new double[size];
	
		//first comma
		pos = geometry.find(cm,lLength);
		size_t pos2 = geometry.find(dlm,lLength);

		//get first point
		lon = geometry.substr(lLength,pos-lLength);
		lat = geometry.substr(pos+1,pos2-pos-1);
		
		stringstream s1(stringstream::in|stringstream::out),s2(stringstream::in|stringstream::out);
		s1 << lon;
		s1 >> points[0];
		s2 << lat;
		s2 >> points[1];

		//cout << "First point of LINE " << geometry << " has Longitude: " << points[0] << " and Latitute: " << points[1] << endl; 

		//first delimeter
		pos = geometry.find(dlm,lLength);
		
		//second comma
		pos2 = geometry.find(cm,pos+1);
		unsigned ind = 2;
		size_t temp;	
	
		//get remaining points
		while(geometry.find(dlm,pos2+1)!=std::string::npos){

			lon = geometry.substr(pos+1,pos2-pos-1);
			pos = geometry.find(dlm,pos2+1);
			lat = geometry.substr (pos2+1,pos-pos2-1);

			stringstream ss1(stringstream::in|stringstream::out),ss2(stringstream::in|stringstream::out);
			ss1 << lon;
			ss1 >> points[ind];
			ind++;
			ss2 << lat;
			ss2 >> points[ind];
			ind++;		
			
			//next comma
			pos2 = geometry.find(cm,pos+1);	
		}

		//get last point
		pos2 = geometry.find(cm,pos+1);
		
		lon = geometry.substr(pos+1,pos2-pos+1);
		lat = geometry.substr (pos2+1,geometry.find(endS,pos2+1)-pos2-1);

		stringstream sss1(stringstream::in|stringstream::out),sss2(stringstream::in|stringstream::out);
		sss1 << lon;
		sss1 >> points[ind];
		ind++;
		sss2 << lat;
		sss2 >> points[ind];
		ind++;

		//cout << "Last point of LINE " << geometry << " has Longitude: " << points[ind-2] << " and Latitute: " << points[ind-1] << endl;

		if(ind!=size){

			cout << "Something is wrong with the line index." << endl;
			cout << "Index: " << ind << " Size: " << size << endl;
			exit(-1);
		}	

		//getchar();

		return 1;
	}
	else if(pl){//polygon, many points
	
		size = 2*countPoints(geometry,dlm); 
		points = new double[size];
	
		//first comma
		pos = geometry.find(cm,plLength);
		size_t pos2 = geometry.find(dlm,plLength);

		//get first point
		lon = geometry.substr(plLength,pos-plLength);
		lat = geometry.substr(pos+1,pos2-pos-1);
		
		stringstream s1(stringstream::in|stringstream::out),s2(stringstream::in|stringstream::out);
		s1 << lon;
		s1 >> points[0];
		s2 << lat;
		s2 >> points[1];

		//cout << "First point of POLYGON " << geometry << " has Longitude: " << points[0] << " and Latitute: " << points[1] << endl;

		//first delimeter
		pos = geometry.find(dlm,plLength);
		
		//second comma
		pos2 = geometry.find(cm,pos+1); 
		unsigned ind = 2;
		size_t temp;	

		//get remaining points
		while(geometry.find(dlm,pos2+1)!=std::string::npos){

			lon = geometry.substr(pos+1,pos2-pos-1);

			//next delimeter
			pos = geometry.find(dlm,pos2+1);

			lat = geometry.substr (pos2+1,pos-pos2-1);

			//cout << "Lon: " << lon << " Lat: " << lat << endl;

			stringstream ss1(stringstream::in|stringstream::out),ss2(stringstream::in|stringstream::out);
			ss1 << lon;
			ss1 >> points[ind];
			ind++;
			ss2 << lat;
			ss2 >> points[ind];
			ind++;		
			
			//next comma
			pos2 = geometry.find(cm,pos+1);	
		}

		//cout << "Index: " << ind << endl;

		//get last point
		pos2 = geometry.find(cm,pos+1);
		
		lon = geometry.substr(pos+1,pos2-pos-1);
		lat = geometry.substr (pos2+1,geometry.find(endS,pos2+1)-pos2-1);

		stringstream sss1(stringstream::in|stringstream::out),sss2(stringstream::in|stringstream::out);
		sss1 << lon;
		sss1 >> points[ind];
		sss2 << lat;
		ind++;
		sss2 >> points[ind];
		ind++;

		//cout << "Last point of POLYGON " << geometry << " has Longitude: " << points[ind-2] << " and Latitute: " << points[ind-1] << endl;

		if(ind!=size){

			cout << "Something is wrong with the polygon index." << endl;
			cout << "Index: " << ind << " Size: " << size << endl;
			exit(-1);
		}

		return 2;	
	}
	else{//multipoint, many points
	
		size = 2*countPoints(geometry,dlm); 
		points = new double[size];
	
		//first comma
		pos = geometry.find(cm,mpLength);
		size_t pos2 = geometry.find(dlm,mpLength);

		//get first point
		lon = geometry.substr(mpLength,pos-mpLength);
		lat = geometry.substr(pos+1,pos2-pos-1);
		
		stringstream s1(stringstream::in|stringstream::out),s2(stringstream::in|stringstream::out);
		s1 << lon;
		s1 >> points[0];
		s2 << lat;
		s2 >> points[1];

		//cout << "First point of MULTIPOINT " << geometry << " has Longitude: " << points[0] << " and Latitute: " << points[1] << endl;

		//first delimeter
		pos = geometry.find(dlm,mpLength);
		
		//second comma
		pos2 = geometry.find(cm,pos+1); 
		unsigned ind = 2;
		size_t temp;	

		//get remaining points
		while(geometry.find(dlm,pos2+1)!=std::string::npos){

			lon = geometry.substr(pos+1,pos2-pos-1);

			//next delimeter
			pos = geometry.find(dlm,pos2+1);

			lat = geometry.substr (pos2+1,pos-pos2-1);

			//cout << "Lon: " << lon << " Lat: " << lat << endl;

			stringstream ss1(stringstream::in|stringstream::out),ss2(stringstream::in|stringstream::out);
			ss1 << lon;
			ss1 >> points[ind];
			ind++;
			ss2 << lat;
			ss2 >> points[ind];
			ind++;		
			
			//next comma
			pos2 = geometry.find(cm,pos+1);	
		}

		//cout << "Index: " << ind << endl;

		//get last point
		pos2 = geometry.find(cm,pos+1);
		
		lon = geometry.substr(pos+1,pos2-pos-1);
		lat = geometry.substr (pos2+1,geometry.find(endS,pos2+1)-pos2-1);

		stringstream sss1(stringstream::in|stringstream::out),sss2(stringstream::in|stringstream::out);
		sss1 << lon;
		sss1 >> points[ind];
		sss2 << lat;
		ind++;
		sss2 >> points[ind];
		ind++;

		//cout << "Last point of MULTIPOINT " << geometry << " has Longitude: " << points[ind-2] << " and Latitute: " << points[ind-1] << endl;

		if(ind!=size){

			cout << "Something is wrong with the multipoint index." << endl;
			cout << "Index: " << ind << " Size: " << size << endl;
			exit(-1);
		}

		return 3;	
	}
*/
}
//---------------------------------------------------------------------------
void Selection::Distance::eval(Result& result)
   // Evaluate the predicate
{	
   // get the ids
   Result r1,r2;
   left->eval(r1);
   right->eval(r2);

   cnt++;

   //bool flag = false;

   if (selection->input->verified) 
   {
		//std::cout << "Distance result is verified..." << std::endl;
		result.setBoolean(true);
		//cout << "left1: " << r1.id << " right1: " << r2.id << endl;	

		//cout << "Vers: " << selection->verif << endl;
		//getchar();

		//flag = true;

		#if SEMIJOIN
			//keep verified inputs in a set
			if (selection->verifiedEntries.insert(r1.id).second) {
			
				selection->verif++;		

				//remove from non-verified if previously inserted
				if ((selection->vIter=selection->nonVerifiedEntries.find(r1.id))!=selection->nonVerifiedEntries.end()) {
			
					selection->vIterStop = selection->nonVerifiedEntries.upper_bound(selection->vIter->first);		
					selection->nonVerifiedEntries.erase(selection->vIter,selection->vIterStop);
				}
			}
		#else
			selection->verif++;
		#endif

		return;
   }

   #if SEMIJOIN
   	//entry is not verified

		if (selection->verifiedEntries.count(r1.id)){//if entry exists just skip it

			//remove from non-verified if previously inserted
			if ((selection->vIter=selection->nonVerifiedEntries.find(r1.id))!=selection->nonVerifiedEntries.end()){
			
				selection->vIterStop = selection->nonVerifiedEntries.upper_bound(selection->vIter->first);		
				selection->nonVerifiedEntries.erase(selection->vIter,selection->vIterStop);
			}

			//continue with the next entry
			result.setBoolean(false);

			return;
		}
		else {  //keep non-verified pairs
			//insert new pair
			selection->nonVerifiedEntries.insert(pair<unsigned,unsigned> (r1.id,r2.id));
			result.setBoolean(false);
			return;
		}
   #endif

   // retrieve geometries
   r1.ensureString(selection);
   r2.ensureString(selection);

   //std::cout << "Entered Distance selection..." << std::endl;
   //cout << "left coords: " << r1.value << " right coords: " << r2.value << endl;
 
   // get coordinates. type of geometry -> 0: point, 1:line, 2:polygon, 3:multipoint
   unsigned type1 = getLongLat(r1.value, points1, size1);
   unsigned type2 = getLongLat(r2.value, points2, size2);

   //cout << "n1: " << n1 << " n2: " << n2 << " n3: " << n3 << " n4: " << n4 << endl;
  

   //cout << "S1: " << size1 << ", S2: " << size2 << endl;
   //cout << "(" << fixed << points1[0] << ", " << fixed << points1[1] << ")" << endl;
   //cout << "(" << fixed << points2[0] << ", " << fixed << points2[1] << ")" << endl;

   bool f = false;

   if (!type1 && !type2) {  // both are points

	   //cout << "Point vs Point" << endl;
	   //cout << "DIST: " << distThreshold << endl;

	   if (type==QueryGraph::Filter::Less) {
	   	//cout << "LESS" << endl;
			//result.setBoolean(((n1-n3)*(n1-n3) + (n2-n4)*(n2-n4))<distThreshold);
			//cout << "DIST: " << distance2(points1,points2) << endl;
			result.setBoolean(distance2(points1,points2)<distThreshold);
	   }
	   else if (type==QueryGraph::Filter::LessOrEqual)
			//result.setBoolean(((n1-n3)*(n1-n3) + (n2-n4)*(n2-n4))<=distThreshold);
			result.setBoolean(distance2(points1,points2)<=distThreshold);
	   else if (type==QueryGraph::Filter::Greater)
			//result.setBoolean(((n1-n3)*(n1-n3) + (n2-n4)*(n2-n4))>distThreshold);
			result.setBoolean(distance2(points1,points2)>distThreshold);
	   else if (type==QueryGraph::Filter::GreaterOrEqual)
			//result.setBoolean(((n1-n3)*(n1-n3) + (n2-n4)*(n2-n4))>=distThreshold);
			result.setBoolean(distance2(points1,points2)>=distThreshold);
	   else if (type==QueryGraph::Filter::Equal)
			//result.setBoolean(((n1-n3)*(n1-n3) + (n2-n4)*(n2-n4))==distThreshold);
			result.setBoolean(distance2(points1,points2)>=distThreshold);
   }
   else if ((f=(!type1 && type2==1)) || (!type2 && type1==1)) {  // point line

		//cout << "Point vs Line" << endl;

		if (f) {  // point line

			result.setBoolean(pointLineStringDist2(points1,points2,size2)<distThreshold);

			//deallocate arrays
			//delete[] points1;
			//delete[] points2;

			return;
		}
		// line point
		result.setBoolean(pointLineStringDist2(points2,points1,size1)<distThreshold);	
   }
   else if ((f=(!type1 && type2==2)) || (!type2 && type1==2)) {  // point polygon

		//cout << "Point vs Polygon" << endl;

		if (f) {  // point polygon

			/*
			double dst = pointPolygonDist2(points1,points2,size2);

			if(dst<distThreshold){
				printf("dst: %.100f\n",dst);
				cout << r1.value << endl;		
				cout << r2.value << endl;
				double minX,maxX,minY,maxY;
				getRegion(points2,size2,minX,maxX,minY,maxY);

			}
			*/

			result.setBoolean(pointPolygonDist2(points1,points2,size2)<distThreshold);

			//deallocate arrays
			//delete[] points1;
			//delete[] points2; 

			return;
		}

		/*
		double dst = pointPolygonDist2(points2,points1,size1);

		if(dst<distThreshold){
			printf("dst1: %.100f\n",dst);
			cout << r1.value << endl;		
			cout << r2.value << endl;

		}
		*/

		// polygon point
		result.setBoolean(pointPolygonDist2(points2,points1,size1)<distThreshold);
   }
   else if ((f=(type1==1 && type2==2)) || (type2==1 && type1==2)) {  // line polygon

		//cout << "Line vs polygon" << endl;

		if (f) {  // line polygon

			result.setBoolean(lineStringPolygonDist2(points1,size1,points2,size2)<distThreshold);

			//deallocate arrays
			//delete[] points1;
			//delete[] points2; 

			return;
		}

		// polygon line
		result.setBoolean(lineStringPolygonDist2(points2,size2,points1,size1)<distThreshold);
   }
   else if (type1==1 && type2==1) {  // line line

		//cout << "Line vs Line" << endl;
		result.setBoolean(lineStringLineStringDist2(points1,size1,points2,size2)<distThreshold);
   }
   else if (type1==2 && type2==2) {  // polygon polygon

		//cout << "Polygon vs Polygon" << endl;
		result.setBoolean(polygonPolygonDist2(points1,points2,size1,size2)<distThreshold);
   }
   else if ((f=(!type1 && type2==3)) || (!type2 && type1==3)) {  // point multipoint

		if (f) {  // point multipoint

			result.setBoolean(pointMultipointDist2(points1,points2,size2)<distThreshold);

			//deallocate arrays
			//delete[] points1;
			//delete[] points2; 

			return;
		}

		// multipoint point
		result.setBoolean(pointMultipointDist2(points2,points1,size1)<distThreshold);
   }
   else if (type1==3 && type2==3) {  // multipoint multipoint

		result.setBoolean(multipointMultipointDist2(points1,points2,size1,size2)<distThreshold);
   }
   else if ((f=(type1==1 && type2==3)) || (type2==1 && type1==3)) {  // line multipoint

		if (f) {  // line multipoint

			result.setBoolean(lineStringMultipointDist2(points1,size1,points2,size2)<distThreshold);

			//deallocate arrays
			//delete[] points1;
			//delete[] points2; 

			return;
		}

		// multipoint line
		result.setBoolean(lineStringMultipointDist2(points2,size2,points1,size1)<distThreshold);
   }
   else if ((f=(type1==2 && type2==3)) || (type2==2 && type1==3)) {  // polygon multipoint

		if (f) {  //polygon multipoint

			result.setBoolean(polygonMultipointDist2(points1,points2,size1,size2)<distThreshold);

			//deallocate arrays
			//delete[] points1;
			//delete[] points2; 

			return;
		}

		// multipoint polygon
		result.setBoolean(polygonMultipointDist2(points2,points1,size2,size1)<distThreshold);
   }
   else {
		cout << "Selection: Could not recognize pair of geometries: " << type1 << " " << type2 << endl;
		exit(-1);
   }

	//deallocate arrays
	//delete[] points1;
	//delete[] points2; 
}
//---------------------------------------------------------------------------
void Selection::Distance::getRegion(double *points, unsigned size, double& minX, double& maxX, double& minY, double& maxY)
   // Get the coordinates of the MBR for a given geoemetry
{
	minX = minY = 1000000;	//just a big value
	maxX = maxY = 0;	//the smallest possible given that the coordinates are normalized in [0,360]

	//cout << "Size: " << size << endl;

	for(unsigned i=0;i<size;i++){

		if(i%2){//latitude

			if(minY>points[i])
				minY = points[i];

			if(maxY<points[i])
				maxY = points[i];
		}
		else{//longitude

			if(minX>points[i])
				minX = points[i];

			if(maxX<points[i])
				maxX = points[i];	
		}
	}

	//printf("minX maxX minY maxY: %.6f %.6f %.6f %.6f\n",minX,maxX,minY,maxY);
	//cout << "minX: " << minX << " maxX: " << maxX << " minY: " << minY << " maxY: " << maxY << endl;
	//getchar();
}
//---------------------------------------------------------------------------
string Selection::Distance::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   stringstream ss(stringstream::in|stringstream::out);
   ss << distThreshold;
   std::string e;
   ss >> e;
   std::string rel;

   if(type==QueryGraph::Filter::Less)
	rel = "<";
   else if(type==QueryGraph::Filter::LessOrEqual)
	rel = "<=";
   else if(type==QueryGraph::Filter::Greater)
	rel = ">";
   else if(type==QueryGraph::Filter::GreaterOrEqual)
	rel = ">=";
   else if(type==QueryGraph::Filter::Equal)
	rel = "=";

   return "(Distance("+left->print(out)+" "+right->print(out)+"))"+rel+"("+e+")";
}
#endif
//---------------------------------------------------------------------------
string Selection::WithinF::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "Within("+input->print(out)+") IN "+geo->toString();
}
//---------------------------------------------------------------------------
string Selection::DistanceF::print(PlanPrinter &out)
   // Print the predicate (debugging only)
{
   stringstream ss(stringstream::in | stringstream::out);
   ss << distThreshold;
   std::string e;
   ss >> e;
   std::string rel;

   if (type == QueryGraph::Filter::Less)
		rel = "<";
   else if (type == QueryGraph::Filter::LessOrEqual)
		rel = "<=";
   else if (type == QueryGraph::Filter::Greater)
		rel = ">";
   else if (type == QueryGraph::Filter::GreaterOrEqual)
		rel = ">=";
   else if (type == QueryGraph::Filter::Equal)
		rel = "=";

   return "(Distance("+arg1->print(out)+" "+arg2->print(out)+"))"+rel+"("+e+")";
}
//---------------------------------------------------------------------------
void Selection::Less::eval(Result& result)
   // Evaluate the predicate
{
   Result l,r;
   left->eval(l);
   right->eval(r);

   //std::cout << "Entered normal Less..." << std::endl;

   //std::cout << "Values before: " << l.value << " " << r.value << std::endl;
   //std::cout << "ID: " << l.hasId() << " String: " << l.hasString() << std::endl;
   //std::cout << "ID2: " << r.hasId() << " String2: " << r.hasString() << std::endl;

   // XXX implement type based comparisons!
   l.ensureString(selection);
   r.ensureString(selection);

   //std::cout << "Values after: " << l.value << " " << r.value << std::endl;
   //std::cout << "ID: " << l.hasId() << " String: " << l.hasString() << std::endl;
   //std::cout << "ID2: " << r.hasId() << " String2: " << r.hasString() << std::endl;

   result.setBoolean(l.value<r.value);
}
//---------------------------------------------------------------------------
string Selection::Less::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "("+left->print(out)+")<("+right->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::LessOrEqual::eval(Result& result)
   // Evaluate the predicate
{
   Result l,r;
   left->eval(l);
   right->eval(r);

   // XXX implement type based comparisons!
   l.ensureString(selection);
   r.ensureString(selection);
   result.setBoolean(l.value<=r.value);
}
//---------------------------------------------------------------------------
string Selection::LessOrEqual::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "("+left->print(out)+")<=("+right->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::Plus::eval(Result& result)
   // Evaluate the predicate
{
   Result l,r;
   left->eval(l);
   right->eval(r);

   l.ensureString(selection);
   r.ensureString(selection);
   stringstream s;
   s << (atof(l.value.c_str())+atof(r.value.c_str()));
   result.setLiteral(s.str());
}
//---------------------------------------------------------------------------
string Selection::Plus::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "("+left->print(out)+")+("+right->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::Minus::eval(Result& result)
   // Evaluate the predicate
{
   Result l,r;
   left->eval(l);
   right->eval(r);

   l.ensureString(selection);
   r.ensureString(selection);
   stringstream s;
   s << (atof(l.value.c_str())-atof(r.value.c_str()));
   result.setLiteral(s.str());
}
//---------------------------------------------------------------------------
string Selection::Minus::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "("+left->print(out)+")-("+right->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::Mul::eval(Result& result)
   // Evaluate the predicate
{
   Result l,r;
   left->eval(l);
   right->eval(r);

   l.ensureString(selection);
   r.ensureString(selection);
   stringstream s;
   s << (atof(l.value.c_str())*atof(r.value.c_str()));
   result.setLiteral(s.str());
}
//---------------------------------------------------------------------------
string Selection::Mul::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "("+left->print(out)+")*("+right->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::Div::eval(Result& result)
   // Evaluate the predicate
{
   Result l,r;
   left->eval(l);
   right->eval(r);

   l.ensureString(selection);
   r.ensureString(selection);
   stringstream s;
   s << (atof(l.value.c_str())/atof(r.value.c_str()));
   result.setLiteral(s.str());
}
//---------------------------------------------------------------------------
string Selection::Div::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "("+left->print(out)+")/("+right->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::Not::eval(Result& result)
   // Evaluate the predicate
{
   result.setBoolean(!input->check());
}
//---------------------------------------------------------------------------
string Selection::Not::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "!"+input->print(out);
}
//---------------------------------------------------------------------------
void Selection::Neg::eval(Result& result)
   // Evaluate the predicate
{
   Result i;
   input->eval(i);

   i.ensureString(selection);
   stringstream s;
   s << (-atof(i.value.c_str()));
   result.setLiteral(s.str());
}
//---------------------------------------------------------------------------
string Selection::Neg::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "-"+input->print(out);
}
//---------------------------------------------------------------------------
void Selection::Null::eval(Result& result)
   // Evaluate the predicate
{
   //the max unsigned integer
   result.setId(~0u);
}
//---------------------------------------------------------------------------
string Selection::Null::print(PlanPrinter& /*out*/)
   // Print the predicate (debugging only)
{
   return "NULL";
}
//---------------------------------------------------------------------------
void Selection::False::eval(Result& result)
   // Evaluate the predicate
{
   result.setBoolean(false);
}
//---------------------------------------------------------------------------
string Selection::False::print(PlanPrinter& /*out*/)
   // Print the predicate (debugging only)
{
   return "false";
}
//---------------------------------------------------------------------------
void Selection::Variable::eval(Result &result)
   // Evaluate the Variable predicate
{
	//cout << "Evaluate Variable predicate" << endl;
   //std::cout << "Setting the value " << reg->value << " of the variable" << std::endl;
   result.setId(reg->value);
}
//---------------------------------------------------------------------------
void Selection::Variable::setResolved(bool f)
   // Pass the info about whether the var is spatial or not to its register so that we can avoid the dictionary lookup 
{
   reg->valueIsResolved = f;
}
//---------------------------------------------------------------------------
string Selection::Variable::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return out.formatRegister(reg);
}
//---------------------------------------------------------------------------
void Selection::ConstantLiteral::eval(Result& result)
   // Evaluate the predicate
{
   //std::cout << "Setting the id " << id << " of the constant." << std::endl;

   result.setId(id);
}
//---------------------------------------------------------------------------
string Selection::ConstantLiteral::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return out.formatValue(id);
}
//---------------------------------------------------------------------------
void Selection::TemporaryConstantLiteral::eval(Result& result)
   // Evaluate the predicate
{

   //std::cout << "Setting the value " << value << " of the temporary constant." << std::endl;

   result.setLiteral(value);
 
}
//---------------------------------------------------------------------------
string Selection::TemporaryConstantLiteral::getValue()
{
   return value;
}
//---------------------------------------------------------------------------
string Selection::TemporaryConstantLiteral::print(PlanPrinter& /*out*/)
   // Print the predicate (debugging only)
{
   return value;
}
//---------------------------------------------------------------------------
void Selection::ConstantIRI::eval(Result& result)
   // Evaluate the predicate
{
   result.setId(id);
}
//---------------------------------------------------------------------------
string Selection::ConstantIRI::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return out.formatValue(id);
}
//---------------------------------------------------------------------------
void Selection::TemporaryConstantIRI::eval(Result& result)
   // Evaluate the predicate
{
   result.setIRI(value);
}
//---------------------------------------------------------------------------
string Selection::TemporaryConstantIRI::print(PlanPrinter& /*out*/)
   // Print the predicate (debugging only)
{
   return value;
}
//---------------------------------------------------------------------------
void Selection::FunctionCall::setSelection(Selection* s)
   // Set the selection
{
   Predicate::setSelection(s);
   for (vector<Predicate*>::iterator iter=args.begin(),limit=args.end();iter!=limit;++iter)
      (*iter)->setSelection(s);
}
//---------------------------------------------------------------------------
void Selection::FunctionCall::eval(Result& result)
   // Evaluate the predicate
{
   result.setId(~0u); // XXX perform the call
}
//---------------------------------------------------------------------------
string Selection::FunctionCall::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   string result="<"+func+">(";
   for (vector<Predicate*>::iterator iter=args.begin(),limit=args.end();iter!=limit;++iter) {
      if (iter!=args.begin())
         result+=",";
      result+=(*iter)->print(out);
   }
   result+=")";
   return result;
}
//---------------------------------------------------------------------------
void Selection::BuiltinStr::eval(Result& result)
   // Evaluate the predicate
{
   input->eval(result);
   result.ensureString(selection);
   result.flags|=Result::typeAvailable;
   result.type=Type::Literal;
}
//---------------------------------------------------------------------------
string Selection::BuiltinStr::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "str("+input->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::BuiltinLang::eval(Result& result)
   // Evaluate the predicate
{
   Result i;
   input->eval(result);
   i.ensureType(selection);

   if (i.type!=Type::CustomLanguage) {
      result.setLiteral("");
   } else {
      result.ensureSubType(selection);
      result.setLiteral(i.subTypeValue);
   }
}
//---------------------------------------------------------------------------
string Selection::BuiltinLang::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "lang("+input->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::BuiltinLangMatches::eval(Result& result)
   // Evaluate the predicate
{
   Result l,r;
   left->eval(l);
   l.ensureType(selection);

   if (l.type!=Type::CustomLanguage) {
      result.setBoolean(false);
      return;
   }
   l.ensureSubType(selection);

   right->eval(r);
   r.ensureString(selection);

   result.setBoolean(l.subTypeValue==r.value); // XXX implement language range checks
}
//---------------------------------------------------------------------------
string Selection::BuiltinLangMatches::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return  "langMatches("+left->print(out)+","+right->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::BuiltinDatatype::eval(Result& result)
   // Evaluate the predicate
{
   Result i;
   input->eval(result);
   i.ensureType(selection);

   switch (i.type) {
      case Type::URI: result.setLiteral("http://www.w3.org/2001/XMLSchema#URI"); break;
      case Type::Literal:
      case Type::CustomLanguage:
      case Type::String: result.setLiteral("http://www.w3.org/2001/XMLSchema#string"); break;
      case Type::Integer: result.setLiteral("http://www.w3.org/2001/XMLSchema#integer"); break;
      case Type::Decimal: result.setLiteral("http://www.w3.org/2001/XMLSchema#decimal"); break;
      case Type::Double: result.setLiteral("http://www.w3.org/2001/XMLSchema#double"); break;
      case Type::Boolean: result.setLiteral("http://www.w3.org/2001/XMLSchema#boolean"); break;
      case Type::CustomType: i.ensureSubType(selection); result.setLiteral(i.subTypeValue); break;
   }
}
//---------------------------------------------------------------------------
string Selection::BuiltinDatatype::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "datatype("+input->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::BuiltinBound::eval(Result& result)
   // Evaluate the predicate
{
   result.setBoolean(~reg->value);
}
//---------------------------------------------------------------------------
string Selection::BuiltinBound::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "bound("+out.formatRegister(reg)+")";
}
//---------------------------------------------------------------------------
void Selection::BuiltinSameTerm::eval(Result& result)
   // Evaluate the predicate
{
   Result l,r;
   left->eval(l);
   right->eval(r);

   // Cheap case
   if (l.hasId()&&r.hasId()) {
      result.setBoolean(l.id==r.id);
      return;
   }

   // Expensive tests
   l.ensureType(selection);
   r.ensureType(selection);
   if ((l.type!=r.type)||(Type::hasSubType(l.type)&&(l.subType!=r.subType))) {
      result.setBoolean(false);
      return;
   }
   l.ensureString(selection);
   r.ensureString(selection);
   result.setBoolean(l.value==r.value);
}
//---------------------------------------------------------------------------
string Selection::BuiltinSameTerm::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "sameTerm("+left->print(out)+","+right->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::BuiltinIsIRI::eval(Result& result)
   // Evaluate the predicate
{
   Result i;
   input->eval(i);
   i.ensureType(selection);
   result.setBoolean(i.type==Type::URI);
}
//---------------------------------------------------------------------------
string Selection::BuiltinIsIRI::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "isIRI("+input->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::BuiltinIsBlank::eval(Result& result)
   // Evaluate the predicate
{
   Result i;
   input->eval(i);
   i.ensureType(selection);
   if (i.type!=Type::URI) {
      result.setBoolean(false);
   } else {
      i.ensureString(selection);
      result.setBoolean(i.value.substr(0,2)=="_:");
   }
}
//---------------------------------------------------------------------------
string Selection::BuiltinIsBlank::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "isBlank("+input->print(out)+")";
}
//---------------------------------------------------------------------------
void Selection::BuiltinIsLiteral::eval(Result& result)
   // Evaluate the predicate
{
   Result i;
   input->eval(i);
   i.ensureType(selection);
   result.setBoolean(i.type!=Type::URI);
}
//---------------------------------------------------------------------------
string Selection::BuiltinIsLiteral::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   return "isLiteral("+input->print(out)+")";
}
//---------------------------------------------------------------------------
Selection::BuiltinRegEx::~BuiltinRegEx()
   // Destructor
{
   delete arg1;
   delete arg2;
   delete arg3;
}
//---------------------------------------------------------------------------
void Selection::BuiltinRegEx::setSelection(Selection* s)
   // Set the selection
{
   Predicate::setSelection(s);
   arg1->setSelection(s);
   arg2->setSelection(s);
   if (arg3) arg3->setSelection(s);
}
//---------------------------------------------------------------------------
void Selection::BuiltinRegEx::eval(Result& result)
   // Evaluate the predicate
{
/*#ifdef CONFIG_TR1
   Result text,pattern;
   arg1->eval(text);
   arg2->eval(pattern);

   try {
      pattern.ensureString(selection);
      std::tr1::regex r(pattern.value.c_str());
      text.ensureString(selection);
      result.setBoolean(std::tr1::regex_match(text.value.begin(),text.value.end(),r));
      return;
   } catch (const std::tr1::regex_error&) {
      result.setBoolean(false);
   }
#else */
   result.setBoolean(false);
//#endif
}
//---------------------------------------------------------------------------
string Selection::BuiltinRegEx::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   string result="regex("+arg1->print(out)+","+arg2->print(out);
   if (arg3) {
      result+=",";
      result+=arg3->print(out);
   }
   result+=")";
   return result;
}
//---------------------------------------------------------------------------
void Selection::BuiltinIn::setSelection(Selection* s)
   // Set the selection
{
   Predicate::setSelection(s);
   probe->setSelection(s);
   for (vector<Predicate*>::iterator iter=args.begin(),limit=args.end();iter!=limit;++iter)
      (*iter)->setSelection(s);
}
//---------------------------------------------------------------------------
void Selection::BuiltinIn::eval(Result& result)
   // Evaluate the predicate
{
   Result p,c;
   probe->eval(p);

   for (vector<Predicate*>::iterator iter=args.begin(),limit=args.end();iter!=limit;++iter) {
      (*iter)->eval(c);
      if (p.hasId()&&c.hasId()) {
         if (p.id==c.id) {
            result.setBoolean(true);
            return;
         }
      } else {
         p.ensureString(selection);
         c.ensureString(selection);
         if (p.value==c.value) {
            result.setBoolean(true);
            return;
         }
      }
   }
   result.setBoolean(false);
}
//---------------------------------------------------------------------------
string Selection::BuiltinIn::print(PlanPrinter& out)
   // Print the predicate (debugging only)
{
   string result="in("+probe->print(out);
   for (vector<Predicate*>::iterator iter=args.begin(),limit=args.end();iter!=limit;++iter) {
      result+=",";
      result+=(*iter)->print(out);
   }
   result+=")";
   return result;
}
//---------------------------------------------------------------------------
Selection::Selection(Operator *input, Runtime &runtime, Predicate *predicate, double expectedOutputCardinality)
   : Operator(expectedOutputCardinality), input(input), runtime(runtime), predicate(predicate)
   // Constructor
{
	#if SEMIJOIN
		done = false;
	#endif
}
//---------------------------------------------------------------------------
Selection::Selection(Operator *entry_RW, vector<Register*> &entryTail_RW, Runtime &runtime, 
						   Geo::Geometry *g, double expectedOutputCardinality)
	: Operator(expectedOutputCardinality), entry_RW(entry_RW), entryTail_RW(entryTail_RW), runtime(runtime)	
	// Constructor to check false positives when R-tree is used to evaluate the WITHIN predicate
{
	if (g->type == 1) {  // Geometry: POINT
		//cout << "POINT found in R-tree WITHIN" << endl;
		geo_RW = new Geo::Point((Geo::Point*) g);
	}	
	else if (g->type == 2) {  // Geometry: RECTANGLE
		//cout << "RECTANGLE found in R-tree WITHIN" << endl;
		geo_RW = new Geo::Rectangle((Geo::Rectangle*) g);
	}	
	geo_RW->normalize();

	#if RTREE_USED_IN_PLAN
		IS = new RetrieveGeoIDsIndexScan(runtime.getDatabase(), Database::Order_Subject_Predicate_Object);
	#endif
}
//---------------------------------------------------------------------------
#if K_NEAREST_NEIGHBOR
Selection::Selection(Operator *entry, vector<Register*> &entryTail, Runtime &runtime, Geo::Geometry *g, 
                     uint32_t k, double expectedOutputCardinality, Grid *gr, bool encodeFlag)
   : Operator(expectedOutputCardinality), entry(entry), entryTail(entryTail), 
	  runtime(runtime), k_nearest_neighbor(k), grid(gr), kNN_encode(encodeFlag) 
   // Constructor kNN
{
	if (g->type == 1) {  // Geometry: POINT
		//cout << "POINT found in kNN selection" << endl;
		geo = new Geo::Point((Geo::Point*) g);
	}
	else if (g->type == 2) {  // Geometry: RECTANGLE
		//cout << "RECTANGLE found in kNN selection" << endl;
		geo = new Geo::Rectangle((Geo::Rectangle*) g);
	}

	if (encodeFlag) { 
		geo->toUnitSquare(grid->minX, grid->minY, grid->maxX, grid->maxY);
		lastDist = 1000000.0;
		cell_cnt = 0;
		gridCells = (grid->CELLS_PER_DIMENSION * grid->CELLS_PER_DIMENSION); 
	}		
	else { 			
		geo->normalize();
	}

	cnt = 0;

	#if RTREE_USED_IN_PLAN
		IS = new RetrieveGeoIDsIndexScan(runtime.getDatabase(), Database::Order_Subject_Predicate_Object);
	#endif
}
#endif
//---------------------------------------------------------------------------
Selection::~Selection()
   // Destructor
{
	#if K_NEAREST_NEIGHBOR
		delete entry;
		delete geo;
		#if RTREE_USED_IN_PLAN
			delete IS;
		#endif	
	#else
		#if RTREE_USED_IN_PLAN
			delete entry_RW;
			delete geo_RW;
			delete IS;
		#else
   		delete predicate;
   		delete input;
   	#endif	
	#endif
}
//---------------------------------------------------------------------------
unsigned Selection::first()
   // Produce the first tuple
{
#if K_NEAREST_NEIGHBOR 

	if (kNN_encode) {

		#if SORTED_ID

			//cout << "first...kNN ENCODING SORTED Selection Evaluation" << endl;
		
			observedOutputCardinality = 0;				
			verif = filtered = 0;	

			// Get the first tuple
			unsigned count = entry->first();
			if (!count) return 0;
	
			// Get the 's' value of triple 's <hasGeometry> geo'	
			unsigned subValue;
			for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {		
				if ( (*iter)->subVar_exist ) {
					subValue = (*iter)->value;  // get and store the subject id
					break;
				}
			}

			
			// DEFINE 'QUERY' RELATED DATA < START >			
			type1 = geo->type;
			//cout << "type1: " << type1 << endl;		
		
			if (type1 == 1) {  // POINT coordinates in the Unit Square
				rx1 = rx2 = ((Geo::Point*) geo)->getX();	
				ry1 = ry2 = ((Geo::Point*) geo)->getY();	
				// Coords Conversion from Unit Square to Normalized	
				nx1 = nx2 = (rx1 * (grid->maxX - grid->minX)) + grid->minX;	
				ny1 = ny2 = (ry1 * (grid->maxY - grid->minY)) + grid->minY;		
			}
			else if (type1 == 2) {  // RECTANGLE coordinates in the Unit Square  
   			rx1 = ((Geo::Rectangle*) geo)->getP1()->getX(); 
   			ry1 = ((Geo::Rectangle*) geo)->getP1()->getY(); 
   			rx2 = ((Geo::Rectangle*) geo)->getP2()->getX(); 
   			ry2 = ((Geo::Rectangle*) geo)->getP2()->getY();		
				// Coords Conversion from Unit Square to Normalized	
				nx1 = (rx1 * (grid->maxX - grid->minX)) + grid->minX;
				ny1 = (ry1 * (grid->maxY - grid->minY)) + grid->minY;
				nx2 = (rx2 * (grid->maxX - grid->minX)) + grid->minX;
				ny2 = (ry2 * (grid->maxY - grid->minY)) + grid->minY;
			}

			//cout << "(rx1, nx1): (" << rx1 << ", " << nx1 << ")" << endl;
			//cout << "(ry1, ny1): (" << ry1 << ", " << ny1 << ")" << endl;
			//cout << "(rx2, nx2): (" << rx2 << ", " << nx2 << ")" << endl;
			//cout << "(ry2, ny2): (" << ry2 << ", " << ny2 << ")" << endl;

			// get the cell id prefix "c_id" and the cell level "level"
			c_id = getHilbertId(rx1, rx2, ry1, ry2, grid->cSide, grid->c2i, grid->CELLS_PER_DIMENSION, grid->b1, grid->b2, level);
			// get the cell gpos "gpos" 
			gpos = getGridPos(c_id, level, grid->size, grid->MAX_LEVEL, grid->b1);

			if (!level) { 
				cout << "QUERY belongs to the Whole Grid" << endl;
				q_ll[0] = 0.0;												    
				q_ll[1] = 0.0;													  
				q_ur[0] = 1.0;  		
				q_ur[1] = 1.0;	 
				////////////////
				cout << "UNDEFINED YET ERROR (1)!!" << endl;
				exit(-1);
			}
			else if (level < grid->MAX_LEVEL) {
				cout << "QUERY is at Upper Level" << endl;
				gpos -= (((unsigned) 1) << grid->b1);
				cX1 = grid->offsets[gpos][0];
				cY1 = grid->offsets[gpos][1];
				cX2 = grid->offsets[gpos][2];
				cY2 = grid->offsets[gpos][3];	
				q_ll[0] = cX1*(grid->cSide);        		
		    	q_ll[1] = cY1*(grid->cSide);        
				q_ur[0] = (cX2 + 1)*(grid->cSide);  
				q_ur[1] = (cY2 + 1)*(grid->cSide);  
				////////////////
				cout << "UNDEFINED YET ERROR (2)!!" << endl;
				exit(-1);
			}
			else {
				//cout << "QUERY is at Bottom Level" << endl;
				cX1 = cX2 = grid->i2c[c_id].first;
				cY1 = cY2 = grid->i2c[c_id].second;
				q_ll[0] = cX1*(grid->cSide);        
		    	q_ll[1] = cY1*(grid->cSide);        
				q_ur[0] = (cX2 + 1)*(grid->cSide);  
				q_ur[1] = (cY2 + 1)*(grid->cSide);  

				q_center[0] = floor(q_ll[0] / grid->cSide) + 0.5;
				q_center[1] = floor(q_ll[1] / grid->cSide) + 0.5;
				q_cell = grid->c2i[(unsigned)(q_ll[0] / grid->cSide)][(unsigned)(q_ll[1] / grid->cSide)]; 
			}			

			//cout << "(Q[0], Q[1]): (" << q_center[0] << ", " << q_center[1] << ")" << endl;
			//cout << "Q_cell: " << q_cell << endl << endl; 
			// DEFINE 'QUERY' RELATED DATA < END >


			// INIT DIST < START >
			double p1[2];
			p1[0] = nx1;  
			p1[1] = ny1;  // QUERY coords

			double p2[2];
			p2[0] = nx1;
			p2[1] = (q_ur[1] * (grid->maxY - grid->minY)) + grid->minY;
			initDist[0] = distance2(p1, p2);  // U
			//cout << "initDist[0]: " << initDist[0] << endl;	
			p2[0] = (q_ur[0] * (grid->maxX - grid->minX)) + grid->minX;
			p2[1] = ny1;
			initDist[1] = distance2(p1, p2);  // R
			//cout << "initDist[1]: " << initDist[1] << endl;	  
			p2[0] = nx1;	
			p2[1] = (q_ll[1] * (grid->maxY - grid->minY)) + grid->minY;
			initDist[2] = distance2(p1, p2);  // D
			//cout << "initDist[2]: " << initDist[2] << endl;
			p2[0] = (q_ll[0] * (grid->maxX - grid->minX)) + grid->minX;		
			p2[1] = ny1;  
			initDist[3] = distance2(p1, p2);  // L
			//cout << "initDist[3]: " << initDist[3] << endl;			
			// INIT DIST < END >


			// GET & PROCESS THE FIRST SPATIAL ENTITY FROM PREVIOUS OPERATORS < START >			
			// Overlook the non-spatial entities
			while ( !getCellId(subValue, c_id, gpos, level) ) {  // get the cell "c_id", "gpos" and "level"
				entry->next();  // get the next tuple 
				// Get the 's' value of triple 's <hasGeometry> geo'
				for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
					if ( (*iter)->subVar_exist ) { 
						subValue = (*iter)->value;  // get and store the subject id
						break;
					} 		
				}
			}

			// Get and store the entity id
			priorityEntry.subject = subValue; 	
	
			// Get and store the 'geo' id and all the 'rest variable' ids 
			for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
				if ( (*iter)->subVar_exist )
					continue; 		
				else if ( (*iter)->geoVar_exist )
					priorityEntry.geometry = (*iter)->value;  // get the ?geo id
				else  
					priorityEntry.restArgs.push_back((*iter)->value);
			}
	

			if (!level) {  // ENTITY belongs to the Whole Grid
				priorityEntry.zone = 0;
				priorityEntry.bottomLevel = 0;	
			}
			else if (level < grid->MAX_LEVEL) {  // ENTITY is at Upper Level
				// Take the REAL 'lower left' and 'upper right' points of the 'upper cell'
				gpos -= (((unsigned) 1) << grid->b1);
				cX1 = grid->offsets[gpos][0];
				cY1 = grid->offsets[gpos][1];
				cX2 = grid->offsets[gpos][2];
				cY2 = grid->offsets[gpos][3];	
				e_ll[0] = cX1*(grid->cSide);        		
		   	e_ll[1] = cY1*(grid->cSide);        
				e_ur[0] = (cX2 + 1)*(grid->cSide);  
				e_ur[1] = (cY2 + 1)*(grid->cSide);
  		
				// FIND THE NEAREST TO THE 'QUERY' BOTTOM CELL PREFIXES OF THE UPPER CELL < START >		
				if (q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // INTERSECTION case
					// select the 'query' cell
					cX1 = q_ll[0] / grid->cSide;  cY1 = q_ll[1] / grid->cSide;
				}
				else if (e_ur[0] <= q_ll[0] && e_ll[1] >= q_ur[1]) {  // UP-LEFT case
					// select the 'low-right' cell
					cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = e_ll[1] / grid->cSide;
				} 
				else if (e_ur[0] <= q_ll[0] && e_ur[1] <= q_ll[1]) {  // LOW-LEFT case
					// select the 'up-right' cell
					cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 
				}
				else if (e_ll[0] >= q_ur[0] && e_ll[1] >= q_ur[1]) {  // UP-RIGHT case
					// select the 'low-left' cell
					cX1 = e_ll[0] / grid->cSide;  cY1 = e_ll[1] / grid->cSide;
				}
				else if (e_ll[0] >= q_ur[0] && e_ur[1] <= q_ll[1]) {  // LOW-RIGHT case  
					// select the 'up-left' cell
					cX1 = e_ll[0] / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 			  	
				}
				else if (e_ur[0] <= q_ll[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // LEFT case
					// select the 'direct-right' cell
					cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = q_ll[1] / grid->cSide;      
				}
				else if (e_ll[0] >= q_ur[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // RIGHT case
					// select the 'direct-left' cell
					cX1 = e_ll[0] / grid->cSide;  cY1 = q_ll[1] / grid->cSide; 				
				}
				else if (e_ll[1] >= q_ur[1] && q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0]) {  // UP case
					// select the 'direct-low' cell
					cX1 = q_ll[0] / grid->cSide;  cY1 = e_ll[1] / grid->cSide;  
				}
				else if (e_ur[1] <= q_ll[1] && q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0]) {  // LOW case
					// select the 'direct-up' cell	
					cX1 = q_ll[0] / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 
				}
				else { 
					cout << "NEAREST BOTTOM CELL ERROR (1)!!" << endl;
					exit(-1);
				}	
				// FIND THE NEAREST TO THE 'QUERY' BOTTOM CELL PREFIXES OF THE UPPER CELL < END >
 
				// Take the ADJUSTED 'lower left' and 'upper right' points  
				// of the nearest to query 'bottom cell' of the 'upper cell'	
				e_ll[0] = cX1 * grid->cSide;
				e_ll[1] = cY1 * grid->cSide;
				e_ur[0] = (cX1 + 1) * grid->cSide;
				e_ur[1] = (cY1 + 1) * grid->cSide;
				// Consider this 'bottom cell' for the kNN evaluation 
				priorityEntry.bottomLevel = 1;
			}
			else {  // ENTITY is at Bottom Level
				cX1 = cX2 = grid->i2c[c_id].first;
				cY1 = cY2 = grid->i2c[c_id].second;
				e_ll[0] = cX1*(grid->cSide);        		
	    		e_ll[1] = cY1*(grid->cSide);        
				e_ur[0] = (cX2 + 1)*(grid->cSide);  
				e_ur[1] = (cY2 + 1)*(grid->cSide);
				priorityEntry.bottomLevel = 1;
			}				


			// Compute the zone
			if (priorityEntry.bottomLevel == 1) {  // ENTITY Is or Adjusted to Bottom Level
				e_center[0] = cX1 + 0.5;  // cell center (x) in terms of cell side			
				e_center[1] = cY1 + 0.5;  // cell center (y) in terms of cell side
				priorityEntry.zone = max(fabs(e_center[0] - q_center[0]), fabs(e_center[1] - q_center[1]));
				//cout << "zone: " << priorityEntry.zone << endl;		
			} 

			// COMPUTE THE 'MINDIST' OF RETRIEVED ENTITY < START >
			if (priorityEntry.zone == 0) {  // 'ZERO ZONE' case			
				priorityEntry.minDist = 0.0;
			}
			else {  // 'NON-ZERO ZONE' case
				// Check all the possible cases	
				if (e_ur[0] <= q_ll[0] && e_ur[1] <= q_ll[1]) {  // LOW-LEFT case
					// select the 'up-right' point (e_ur[0], e_ur[1])  
					p2[0] = (e_ur[0] * (grid->maxX - grid->minX)) + grid->minX; 	
					p2[1] = (e_ur[1] * (grid->maxY - grid->minY)) + grid->minY;
					priorityEntry.minDist = distance2(p1, p2);
				}
				else if (e_ll[0] >= q_ur[0] && e_ur[1] <= q_ll[1]) {  // LOW-RIGHT case
					// select the 'up-left' point (e_ll[0], e_ur[1])
					p2[0] = (e_ll[0] * (grid->maxX - grid->minX)) + grid->minX;	
					p2[1] = (e_ur[1] * (grid->maxY - grid->minY)) + grid->minY;
					priorityEntry.minDist = distance2(p1, p2);		
				}  	
				else if (e_ur[0] <= q_ll[0] && e_ll[1] >= q_ur[1]) {  // UP-LEFT case
					// select the 'low-right' point (e_ur[0], e_ll[1])
					p2[0] = (e_ur[0] * (grid->maxX - grid->minX)) + grid->minX; 		
					p2[1] = (e_ll[1] * (grid->maxY - grid->minY)) + grid->minY;
					priorityEntry.minDist = distance2(p1, p2); 
				}
				else if (e_ll[0] >= q_ur[0] && e_ll[1] >= q_ur[1]) {  // UP-RIGHT case
					// select the 'low-left' point (e_ll[0], e_ll[1])
					p2[0] = (e_ll[0] * (grid->maxX - grid->minX)) + grid->minX;
					p2[1] = (e_ll[1] * (grid->maxY - grid->minY)) + grid->minY;	
					priorityEntry.minDist = distance2(p1, p2);		
				}
				else if (e_ll[0] == q_ll[0] && e_ll[1] > q_ll[1]) {  // DIRECT-UP case
					p2[0] = nx1;
					p2[1] = q_ur[1] + ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'y' coord 
					p2[1] = (p2[1] * (grid->maxY - grid->minY)) + grid->minY;  // normalized 'y' coord
					priorityEntry.minDist = distance2(p1, p2);	
				}
				else if (e_ll[0] > q_ll[0] && e_ll[1] == q_ll[1]) {  // DIRECT-RIGHT case 
					p2[0] = q_ur[0] + ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'x' coord
					p2[0] = (p2[0] * (grid->maxX - grid->minX)) + grid->minX;  // normalized 'x' coord
					p2[1] = ny1;
					priorityEntry.minDist = distance2(p1, p2);
				}		
				else if (e_ll[0] == q_ll[0] && e_ll[1] < q_ll[1]) {  // DIRECT-DOWN case
					p2[0] = nx1;
					p2[1] = q_ll[1] - ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'y' coord
					p2[1] = (p2[1] * (grid->maxY - grid->minY)) + grid->minY;  // normalized 'y' coord	
					priorityEntry.minDist = distance2(p1, p2);
				}	
				else if (e_ll[0] < q_ll[0] && e_ll[1] == q_ll[1]) {  // DIRECT-LEFT case
					p2[0] = q_ll[0] - ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'x' coord
					p2[0] = (p2[0] * (grid->maxX - grid->minX)) + grid->minX;  // normalized 'x' coord		
					p2[1] = ny1;  
					priorityEntry.minDist = distance2(p1, p2);
				}
				else {
					cout << "MINDIST ERROR (1)!!" << endl;
					exit(-1);
				}
			} 			
			// COMPUTE THE 'MINDIST' OF RETRIEVED ENTITY < END >			
	
			//cout << "minDist: " << priorityEntry.minDist << endl;
			priorityEntry.dir = 4;  // cell	
			H.push(priorityEntry);
			// GET & PROCESS THE FIRST SPATIAL ENTITY FROM PREVIOUS OPERATORS < END >


			unsigned limits[14];  // the 'limit entity conditions' for each level
			unsigned maxLimit;  // the maximum 'limit entity condition' among all
									  // levels (including the Whole Grid level)	
			unsigned prevMaxLimit;  // it is updated only when 'maxLimit' gets higher	
			unsigned maxCell;  // the maximum cell of 'respective zone' (in bottom level)		

			// FIND THE 'LIMIT ENTITY CONDITION' OF WHOLE GRID ('L_0') < START >
			limits[13] = 0 << 31;
			if (grid->grid[grid->size - 1] != 0) {	
				limits[13] |= 1;
				limits[13] |= (0 << 1);
				limits[13] |= ((grid->grid[grid->size - 1] - 1) << (grid->b2 - 1));
			}
			// FIND THE 'LIMIT ENTITY CONDITION' OF WHOLE GRID ('L_0') < END >

			// FIND THE 'LIMIT ENTITY CONDITIONS' OF ZERO ZONE < START >
			maxCell = q_cell;
			findEntityLimits(limits, grid->grid, maxCell, grid->b2);
			// FIND THE 'LIMIT ENTITY CONDITIONS' OF ZERO ZONE < END >

			// FIND THE 'maxLimit' AMONG 'ZERO ZONE' & 'WHOLE GRID' < START >
			maxLimit = limits[0];
			if (limits[1] > maxLimit)   maxLimit = limits[1];
			if (limits[2] > maxLimit)   maxLimit = limits[2];
			if (limits[3] > maxLimit)   maxLimit = limits[3];	
			if (limits[4] > maxLimit)   maxLimit = limits[4];
			if (limits[5] > maxLimit)   maxLimit = limits[5];
			if (limits[6] > maxLimit)   maxLimit = limits[6];
			if (limits[7] > maxLimit)   maxLimit = limits[7];
			if (limits[8] > maxLimit)   maxLimit = limits[8];
			if (limits[9] > maxLimit)   maxLimit = limits[9];
			if (limits[10] > maxLimit)  maxLimit = limits[10];
			if (limits[11] > maxLimit)  maxLimit = limits[11];
			if (limits[12] > maxLimit)  maxLimit = limits[12];
			if (limits[13] > maxLimit)  maxLimit = limits[13];
			prevMaxLimit = maxLimit;	
			// FIND THE 'maxLimit' AMONG 'ZERO ZONE' & 'WHOLE GRID' < END >	
	

			// GET & PROCESS ALL THE SPATIAL ENTITIES FROM PREVIOUS OPERATORS, BASED ON THE WHOLE GRID 
			// & ZERO ZONE 'LIMIT ENTITY CONDITIONS'. THIS IS A PROCEDURE THAT CANNOT BE AVOIDED < START >				
			while ( (count = entry->next()) != 0 ) {
				// Get the 's' value of triple 's <hasGeometry> geo'
				for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
					if ( (*iter)->subVar_exist ) { 
						subValue = (*iter)->value;  // get and store the subject id		
						break;
					} 		
				}

				// Overlook the non-spatial entities
				while ( !getCellId(subValue, c_id, gpos, level) ) {  // get the cell "c_id", "gpos" and "level"
					// Stop taking more entities from previous operators ('non-spatial entity' case)
					if (subValue > maxLimit) {
						count = 0;
						break;
					}
					count = entry->next();  // get the next tuple
					if (!count) break;
					// Get the 's' value of triple 's <hasGeometry> geo'
					for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
						if ( (*iter)->subVar_exist ) { 
							subValue = (*iter)->value;  // get and store the subject id
							break;
						} 		
					}
				}

				if (count != 0) {  // Edit all the spatial entities

					// Clear previous 'argument values'
					priorityEntry.restArgs.clear();

					// Get and store the entity id
					priorityEntry.subject = subValue;
	
					// Get and store the 'geo' id and all the 'rest variable' ids 
					for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
						if ( (*iter)->subVar_exist )
							continue; 		
						else if ( (*iter)->geoVar_exist )
							priorityEntry.geometry = (*iter)->value;  // get the ?geo id
						else 
							priorityEntry.restArgs.push_back((*iter)->value);
					}

					
					if (!level) {  // ENTITY belongs to the Whole Grid  
						priorityEntry.zone = 0; 	
						priorityEntry.bottomLevel = 0; 
					}
					else if (level < grid->MAX_LEVEL) {  // ENTITY is at Upper Level
						// Take the REAL 'lower left' and 'upper right' points of the 'upper cell'
						gpos -= (((unsigned) 1) << grid->b1);
						cX1 = grid->offsets[gpos][0];
						cY1 = grid->offsets[gpos][1];
						cX2 = grid->offsets[gpos][2];
						cY2 = grid->offsets[gpos][3];	
						e_ll[0] = cX1*(grid->cSide);       		
		   			e_ll[1] = cY1*(grid->cSide);        
						e_ur[0] = (cX2 + 1)*(grid->cSide); 
						e_ur[1] = (cY2 + 1)*(grid->cSide); 

						// FIND THE NEAREST TO THE 'QUERY' BOTTOM CELL PREFIXES OF THE UPPER CELL < START >		
						if (q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // INTERSECTION case
							// select the 'query' cell
							cX1 = q_ll[0] / grid->cSide;  cY1 = q_ll[1] / grid->cSide;
						}
						else if (e_ur[0] <= q_ll[0] && e_ll[1] >= q_ur[1]) {  // UP-LEFT case
							// select the 'low-right' cell
							cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = e_ll[1] / grid->cSide;
						} 
						else if (e_ur[0] <= q_ll[0] && e_ur[1] <= q_ll[1]) {  // LOW-LEFT case
							// select the 'up-right' cell
							cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 
						}
						else if (e_ll[0] >= q_ur[0] && e_ll[1] >= q_ur[1]) {  // UP-RIGHT case
							// select the 'low-left' cell
							cX1 = e_ll[0] / grid->cSide;  cY1 = e_ll[1] / grid->cSide;
						}
						else if (e_ll[0] >= q_ur[0] && e_ur[1] <= q_ll[1]) {  // LOW-RIGHT case  
							// select the 'up-left' cell
							cX1 = e_ll[0] / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 			  	
						}
						else if (e_ur[0] <= q_ll[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // LEFT case
							// select the 'direct-right' cell
							cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = q_ll[1] / grid->cSide;      
						}
						else if (e_ll[0] >= q_ur[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // RIGHT case
							// select the 'direct-left' cell
							cX1 = e_ll[0] / grid->cSide;  cY1 = q_ll[1] / grid->cSide; 				
						}
						else if (e_ll[1] >= q_ur[1] && q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0]) {  // UP case
							// select the 'direct-low' cell
							cX1 = q_ll[0] / grid->cSide;  cY1 = e_ll[1] / grid->cSide;  
						}
						else if (e_ur[1] <= q_ll[1] && q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0]) {  // LOW case
							// select the 'direct-up' cell	
							cX1 = q_ll[0] / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 
						}
						else { 
							cout << "NEAREST BOTTOM CELL ERROR (2)!!" << endl;
							exit(-1);
						}	
						// FIND THE NEAREST TO THE 'QUERY' BOTTOM CELL PREFIXES OF THE UPPER CELL < END >

						// Take the ADJUSTED 'lower left' and 'upper right' points  
						// of the nearest to query 'bottom cell' of the 'upper cell'	
						e_ll[0] = cX1 * grid->cSide;
						e_ll[1] = cY1 * grid->cSide;
						e_ur[0] = (cX1 + 1) * grid->cSide;
						e_ur[1] = (cY1 + 1) * grid->cSide;
						// Consider this 'bottom cell' for the kNN evaluation 
						priorityEntry.bottomLevel = 1;
					}
					else {  // ENTITY is at Bottom Level    
						cX1 = cX2 = grid->i2c[c_id].first;
						cY1 = cY2 = grid->i2c[c_id].second;
						e_ll[0] = cX1*(grid->cSide);        		
	    				e_ll[1] = cY1*(grid->cSide);        
						e_ur[0] = (cX2 + 1)*(grid->cSide);  
						e_ur[1] = (cY2 + 1)*(grid->cSide);
						priorityEntry.bottomLevel = 1;
					}								


					// Compute the zone
					if (priorityEntry.bottomLevel == 1) {  // ENTITY Is or Adjusted to Bottom Level
						e_center[0] = cX1 + 0.5;  // cell center (x) in terms of cell side			
						e_center[1] = cY1 + 0.5;  // cell center (y) in terms of cell side
						priorityEntry.zone = max(fabs(e_center[0] - q_center[0]), fabs(e_center[1] - q_center[1]));
						//cout << "zone: " << priorityEntry.zone << endl;	
					}	

					// COMPUTE THE 'MINDIST' OF RETRIEVED ENTITY < START >
					if (priorityEntry.zone == 0) {  // 'ZERO ZONE' case
						priorityEntry.minDist = 0.0;	
					}
					else {  // 'NON-ZERO ZONE' case
						// Check all the possible cases	
						if (e_ur[0] <= q_ll[0] && e_ur[1] <= q_ll[1]) {  // LOW-LEFT case
							// select the 'up-right' point (e_ur[0], e_ur[1])  
							p2[0] = (e_ur[0] * (grid->maxX - grid->minX)) + grid->minX;
							p2[1] = (e_ur[1] * (grid->maxY - grid->minY)) + grid->minY; 
							priorityEntry.minDist = distance2(p1, p2); 
						}
						else if (e_ll[0] >= q_ur[0] && e_ur[1] <= q_ll[1]) {  // LOW-RIGHT case
							// select the 'up-left' point (e_ll[0], e_ur[1])
							p2[0] = (e_ll[0] * (grid->maxX - grid->minX)) + grid->minX;
							p2[1] = (e_ur[1] * (grid->maxY - grid->minY)) + grid->minY;
							priorityEntry.minDist = distance2(p1, p2);	
						}  	
						else if (e_ur[0] <= q_ll[0] && e_ll[1] >= q_ur[1]) {  // UP-LEFT case
							// select the 'low-right' point (e_ur[0], e_ll[1])
							p2[0] = (e_ur[0] * (grid->maxX - grid->minX)) + grid->minX;
							p2[1] = (e_ll[1] * (grid->maxY - grid->minY)) + grid->minY;	
							priorityEntry.minDist = distance2(p1, p2); 	
						}
						else if (e_ll[0] >= q_ur[0] && e_ll[1] >= q_ur[1]) {  // UP-RIGHT case
							// select the 'low-left' point (e_ll[0], e_ll[1])
							p2[0] = (e_ll[0] * (grid->maxX - grid->minX)) + grid->minX;
							p2[1] = (e_ll[1] * (grid->maxY - grid->minY)) + grid->minY;	
							priorityEntry.minDist = distance2(p1, p2);	
						}
						else if (e_ll[0] == q_ll[0] && e_ll[1] > q_ll[1]) {  // DIRECT-UP case
							p2[0] = nx1;
							p2[1] = q_ur[1] + ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'y' coord 
							p2[1] = (p2[1] * (grid->maxY - grid->minY)) + grid->minY;  // normalized 'y' coord
							priorityEntry.minDist = distance2(p1, p2);	
						}
						else if (e_ll[0] > q_ll[0] && e_ll[1] == q_ll[1]) {  // DIRECT-RIGHT case 
							p2[0] = q_ur[0] + ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'x' coord
							p2[0] = (p2[0] * (grid->maxX - grid->minX)) + grid->minX;  // normalized 'x' coord
							p2[1] = ny1;
							priorityEntry.minDist = distance2(p1, p2);
						}	
						else if (e_ll[0] == q_ll[0] && e_ll[1] < q_ll[1]) {  // DIRECT-DOWN case
							p2[0] = nx1;
							p2[1] = q_ll[1] - ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'y' coord
							p2[1] = (p2[1] * (grid->maxY - grid->minY)) + grid->minY;  // normalized 'y' coord	
							priorityEntry.minDist = distance2(p1, p2);
						}	
						else if (e_ll[0] < q_ll[0] && e_ll[1] == q_ll[1]) {  // DIRECT-LEFT case
							p2[0] = q_ll[0] - ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'x' coord
							p2[0] = (p2[0] * (grid->maxX - grid->minX)) + grid->minX;  // normalized 'x' coord		
							p2[1] = ny1;  
							priorityEntry.minDist = distance2(p1, p2);
						}
						else {
							cout << "MINDIST ERROR (2)!!" << endl;
							exit(-1);
						}			
					}
					// COMPUTE THE 'MINDIST' OF RETRIEVED ENTITY < END >

					//cout << "minDist: " << priorityEntry.minDist << endl;
					priorityEntry.dir = 4;  // cell
					H.push(priorityEntry);				
						
					// Stop taking more entities from previous operators ('spatial entity' case)
					if (subValue >= maxLimit) { 
						prevMaxLimit = subValue; 
				   	break;
					}	
				}
				else {  // All the entities are taken OR a non-spatial one is found that is above 'maxLimit'
				   break;
				}
			}		
			// GET & PROCESS ALL THE SPATIAL ENTITIES FROM PREVIOUS OPERATORS, BASED ON THE WHOLE GRID 
			// & ZERO ZONE 'LIMIT ENTITY CONDITIONS'. THIS IS A PROCEDURE THAT CANNOT BE AVOIDED < END >


			// PUT THE 'FIRST ZONE' VIRTUAL RECTANGLES IN 'H' < START >
			priorityEntry.zone = 1;  priorityEntry.dir = 0;  priorityEntry.minDist = initDist[0];  H.push(priorityEntry);  // U
			priorityEntry.zone = 1;  priorityEntry.dir = 1;  priorityEntry.minDist = initDist[1];  H.push(priorityEntry);  // R
			priorityEntry.zone = 1;  priorityEntry.dir = 2;  priorityEntry.minDist = initDist[2];  H.push(priorityEntry);  // D
			priorityEntry.zone = 1;  priorityEntry.dir = 3;  priorityEntry.minDist = initDist[3];  H.push(priorityEntry);  // L
			// PUT THE 'FIRST ZONE' VIRTUAL RECTANGLES IN 'H' < END >


			// QUERY COORDS INITIALIZATION FOR DISTANCE FUNCTIONS < START >
			if (type1 == 1) {  // query is POINT
				size1 = 2;
				points1 = new double[size1];
				points1[0] = nx1;  points1[1] = ny1;
				type1 = 0;  // adjusted for DISTANCE functions
			}	
			else if (type1 == 2) {  // query is RECTANGLE
				size1 = 4;
				points1 = new double[size1];
				points1[0] = nx1;  points1[1] = ny1;
				points1[2] = nx2;  points1[3] = ny2;				
			}
			// QUERY COORDS INITIALIZATION FOR DISTANCE FUNCTIONS < END >
	

			// APPLY CPM < START >
			Result geoRes;
			bool f;

			while (!H.empty()) {

				priorityEntry = H.top();
				if (priorityEntry.minDist >= lastDist)		
					break; 
				// Remove the retrieved Entity from 'H'
				H.pop();  

				if (priorityEntry.dir == 4) {  // 'CELL ENTITY' case	
		
               geoRes.flags = 0;  // reset the flags
					geoRes.flags = Result::idAvailable;		
					geoRes.id = priorityEntry.geometry;  // get the ?geo id
					geoRes.ensureString(this);  // get the ?geo string value ('retrieve geometry')

					type2 = getLongLat(geoRes.value, points2, size2);  // set 'points2' & 'size2' values (NORMALIZED COORDS)
				
					// COMPUTE REAL DISTANCE < START >
					f = false;
		
					if (!type1 && !type2) {  // both are points
						//cout << "Point vs Point" << endl;
						priorityEntry.realDist = distance2(points1, points2);
					}
					else if ((f=(!type1 && type2==1)) || (!type2 && type1==1)) {  // point line
						if (f) {  // point line
							//cout << "Point vs Line" << endl;
							priorityEntry.realDist = pointLineStringDist2(points1, points2, size2);
						}
						else {  // line point
							//cout << "Line vs Point" << endl;
							priorityEntry.realDist = pointLineStringDist2(points2, points1, size1);
						}
					}
   				else if ((f=(!type1 && type2==2)) || (!type2 && type1==2)) {  // point polygon
						if (f) {  // point polygon
							//cout << "Point vs Polygon" << endl;
							priorityEntry.realDist = pointPolygonDist2(points1, points2, size2);
						}
						else {  // polygon point
							//cout << "Polygon vs Point" << endl;
							priorityEntry.realDist = pointPolygonDist2(points2, points1, size1);
						}	
   				}
					else if ((f=(type1==1 && type2==2)) || (type2==1 && type1==2)) {  // line polygon
						if (f) {  // line polygon
							//cout << "Line vs Polygon" << endl;
							priorityEntry.realDist = lineStringPolygonDist2(points1, size1, points2, size2);
						}
						else {  // polygon line
							//cout << "Polygon vs Line" << endl;
							priorityEntry.realDist = lineStringPolygonDist2(points2, size2, points1, size1);
						}
					}	
					else if (type1==1 && type2==1) {  // line line
						//cout << "Line vs Line" << endl;	
						priorityEntry.realDist = lineStringLineStringDist2(points1, size1, points2, size2);
					}
					else if (type1==2 && type2==2) {  // polygon polygon
						//cout << "Polygon vs Polygon" << endl;
						priorityEntry.realDist = polygonPolygonDist2(points1, points2, size1, size2);
					}
					else if ((f=(!type1 && type2==3)) || (!type2 && type1==3)) {  // point multipoint
						if (f) {  // point multipoint	
							//cout << "Point vs Multipoint" << endl;
							priorityEntry.realDist = pointMultipointDist2(points1, points2, size2);
						}
						else {  // multipoint point
							//cout << "Multipoint vs Point" << endl;
							priorityEntry.realDist = pointMultipointDist2(points2, points1, size1);
						}
					}
					else if (type1==3 && type2==3) {  // multipoint multipoint
						//cout << "Multipoint vs Multipoint" << endl;
						priorityEntry.realDist = multipointMultipointDist2(points1, points2, size1, size2);	
					}
					else if ((f=(type1==1 && type2==3)) || (type2==1 && type1==3)) {  // line multipoint
						if (f) {  // line multipoint
							//cout << "Line vs Multipoint" << endl;
							priorityEntry.realDist = lineStringMultipointDist2(points1, size1, points2, size2);
						}
						else {  // multipoint line
							//cout << "Multipoint vs Line" << endl;
							priorityEntry.realDist = lineStringMultipointDist2(points2, size2, points1, size1);
						}
					}
					else if ((f=(type1==2 && type2==3)) || (type2==2 && type1==3)) {  // polygon multipoint
						if (f) {  // polygon multipoint
							//cout << "Polygon vs Multipoint" << endl;
							priorityEntry.realDist = polygonMultipointDist2(points1, points2, size1, size2);
						}
						else {  // multipoint polygon
							//cout << "Multipoint vs Polygon" << endl;
							priorityEntry.realDist = polygonMultipointDist2(points2, points1, size2, size1);	
						}
					}
					else {
						cout << "kNN Selection ERROR: Could not recognize pair of geometries: " << type1 << " " << type2 << endl;
						exit(-1);
					}
					// COMPUTE REAL DISTANCE < END >

					cell_cnt++;
					if (cell_cnt > k_nearest_neighbor) {
						R.push(priorityEntry);				
						R.pop();
					}
					else {  
						R.push(priorityEntry);
					}

					if (cell_cnt >= k_nearest_neighbor) { 
						priorityEntry = R.top();
						lastDist = priorityEntry.realDist;
					}
				}

				else {  // 'RECTANGLE ENTITY' case

					// FIND THE '(priorityEntry.zone, priorityEntry.dir)' 'LIMIT ENTITY CONDITIONS' < START > 
					maxCell = findMaxCellInRecZone(priorityEntry.zone, priorityEntry.dir, grid->c2i, q_ll, grid->cSide); 
					if (maxCell != gridCells)  // '(priorityEntry.zone, priorityEntry.dir)' is WITHIN Grid
						findEntityLimits(limits, grid->grid, maxCell, grid->b2);
					// FIND THE '(priorityEntry.zone, priorityEntry.dir)' 'LIMIT ENTITY CONDITIONS' < END >

					if (maxCell != gridCells) {  // '(priorityEntry.zone, priorityEntry.dir)' is WITHIN Grid

						// PUT THE 'priorityEntry.zone + 1' VIRTUAL RECTANGLES IN 'H' < START >
						priorityEntry.zone = priorityEntry.zone + 1; 
						if (priorityEntry.dir == 0) {  // U
							p2[0] = nx1;
							p2[1] = q_ur[1] + ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'y' coord 
							p2[1] = (p2[1] * (grid->maxY - grid->minY)) + grid->minY;  // normalized 'y' coord
							priorityEntry.minDist = distance2(p1, p2); 
						}	
						else if (priorityEntry.dir == 1) {  // R 
							p2[0] = q_ur[0] + ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'x' coord
							p2[0] = (p2[0] * (grid->maxX - grid->minX)) + grid->minX;  // normalized 'x' coord
							p2[1] = ny1;
							priorityEntry.minDist = distance2(p1, p2);	
						}
						else if (priorityEntry.dir == 2) {  // D	
							p2[0] = nx1;
							p2[1] = q_ll[1] - ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'y' coord
							p2[1] = (p2[1] * (grid->maxY - grid->minY)) + grid->minY;  // normalized 'y' coord	
							priorityEntry.minDist = distance2(p1, p2);
						}
						else if (priorityEntry.dir == 3) {  // L 
							p2[0] = q_ll[0] - ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'x' coord
							p2[0] = (p2[0] * (grid->maxX - grid->minX)) + grid->minX;  // normalized 'x' coord		
							p2[1] = ny1;  
							priorityEntry.minDist = distance2(p1, p2);
						}
						else {
							cout << "DIR ERROR!!" << endl;
							exit(-1);
						}  
						H.push(priorityEntry);
						// PUT THE 'priorityEntry.zone + 1' VIRTUAL RECTANGLES IN 'H' < END >


						// FIND THE 'maxLimit' OF '(priorityEntry.zone, priorityEntry.dir)' < START >
						maxLimit = limits[0];
						if (limits[1] > maxLimit)   maxLimit = limits[1];
						if (limits[2] > maxLimit)   maxLimit = limits[2];
						if (limits[3] > maxLimit)   maxLimit = limits[3];	
						if (limits[4] > maxLimit)   maxLimit = limits[4];
						if (limits[5] > maxLimit)   maxLimit = limits[5];
						if (limits[6] > maxLimit)   maxLimit = limits[6];
						if (limits[7] > maxLimit)   maxLimit = limits[7];
						if (limits[8] > maxLimit)   maxLimit = limits[8];
						if (limits[9] > maxLimit)   maxLimit = limits[9];
						if (limits[10] > maxLimit)  maxLimit = limits[10];
						if (limits[11] > maxLimit)  maxLimit = limits[11];
						if (limits[12] > maxLimit)  maxLimit = limits[12];
						// FIND THE 'maxLimit' OF '(priorityEntry.zone, priorityEntry.dir)' < END >
			

						if (maxLimit > prevMaxLimit)
						{
					
						// GET & PROCESS ALL THE SPATIAL ENTITIES FROM PREVIOUS OPERATORS, BASED ON THE 
						// '(priorityEntry.zone, priorityEntry.dir)' 'LIMIT ENTITY CONDITIONS' < START >				
						while ( (count = entry->next()) != 0 ) {
							// Get the 's' value of triple 's <hasGeometry> geo'
							for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
								if ( (*iter)->subVar_exist ) { 
									subValue = (*iter)->value;  // get and store the subject id		
									break;
								} 		
							}

							// Overlook the non-spatial entities
							while ( !getCellId(subValue, c_id, gpos, level) ) {  // get the cell "c_id", "gpos" and "level"
								// Stop taking more entities from previous operators ('non-spatial entity' case)
								if (subValue > maxLimit) {
									count = 0;
									break;
								}
								count = entry->next();  // get the next tuple
								if (!count) break;
								// Get the 's' value of triple 's <hasGeometry> geo'
								for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
									if ( (*iter)->subVar_exist ) { 
										subValue = (*iter)->value;  // get and store the subject id
										break;
									} 		
								}
							}

							if (count != 0) {  // Edit all the spatial entities

								// Clear previous 'argument values'
								priorityEntry.restArgs.clear();

								// Get and store the entity id
								priorityEntry.subject = subValue;

								// Get and store the 'geo' id and all the 'rest variable' ids 
								for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
									if ( (*iter)->subVar_exist )
										continue; 		
									else if ( (*iter)->geoVar_exist )
										priorityEntry.geometry = (*iter)->value;  // get the ?geo id
									else 
										priorityEntry.restArgs.push_back((*iter)->value);
								}

					
								if (!level) {  // ENTITY belongs to the Whole Grid  
									priorityEntry.zone = 0; 	
									priorityEntry.bottomLevel = 0; 
								}
								else if (level < grid->MAX_LEVEL) {  // ENTITY is at Upper Level
									// Take the REAL 'lower left' and 'upper right' points of the 'upper cell'
									gpos -= (((unsigned) 1) << grid->b1);
									cX1 = grid->offsets[gpos][0];
									cY1 = grid->offsets[gpos][1];
									cX2 = grid->offsets[gpos][2];
									cY2 = grid->offsets[gpos][3];	
									e_ll[0] = cX1*(grid->cSide);       		
		   						e_ll[1] = cY1*(grid->cSide);        
									e_ur[0] = (cX2 + 1)*(grid->cSide); 
									e_ur[1] = (cY2 + 1)*(grid->cSide); 

									// FIND THE NEAREST TO THE 'QUERY' BOTTOM CELL PREFIXES OF THE UPPER CELL < START >		
									if (q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // INTERSECTION case
										// select the 'query' cell
										cX1 = q_ll[0] / grid->cSide;  cY1 = q_ll[1] / grid->cSide;
									}
									else if (e_ur[0] <= q_ll[0] && e_ll[1] >= q_ur[1]) {  // UP-LEFT case
										// select the 'low-right' cell
										cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = e_ll[1] / grid->cSide;
									} 
									else if (e_ur[0] <= q_ll[0] && e_ur[1] <= q_ll[1]) {  // LOW-LEFT case
										// select the 'up-right' cell
										cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 
									}
									else if (e_ll[0] >= q_ur[0] && e_ll[1] >= q_ur[1]) {  // UP-RIGHT case
										// select the 'low-left' cell
										cX1 = e_ll[0] / grid->cSide;  cY1 = e_ll[1] / grid->cSide;
									}
									else if (e_ll[0] >= q_ur[0] && e_ur[1] <= q_ll[1]) {  // LOW-RIGHT case  
										// select the 'up-left' cell
										cX1 = e_ll[0] / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 			  	
									}
									else if (e_ur[0] <= q_ll[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // LEFT case
										// select the 'direct-right' cell
										cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = q_ll[1] / grid->cSide;      
									}
									else if (e_ll[0] >= q_ur[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // RIGHT case
										// select the 'direct-left' cell
										cX1 = e_ll[0] / grid->cSide;  cY1 = q_ll[1] / grid->cSide; 				
									}
									else if (e_ll[1] >= q_ur[1] && q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0]) {  // UP case
										// select the 'direct-low' cell
										cX1 = q_ll[0] / grid->cSide;  cY1 = e_ll[1] / grid->cSide;  
									}
									else if (e_ur[1] <= q_ll[1] && q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0]) {  // LOW case
										// select the 'direct-up' cell	
										cX1 = q_ll[0] / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 
									}
									else { 
										cout << "NEAREST BOTTOM CELL ERROR (3)!!" << endl;
										exit(-1);
									}	
									// FIND THE NEAREST TO THE 'QUERY' BOTTOM CELL PREFIXES OF THE UPPER CELL < END >

									// Take the ADJUSTED 'lower left' and 'upper right' points  
									// of the nearest to query 'bottom cell' of the 'upper cell'	
									e_ll[0] = cX1 * grid->cSide;
									e_ll[1] = cY1 * grid->cSide;
									e_ur[0] = (cX1 + 1) * grid->cSide;
									e_ur[1] = (cY1 + 1) * grid->cSide;
									// Consider this 'bottom cell' for the kNN evaluation 
									priorityEntry.bottomLevel = 1;
								}
								else {  // ENTITY is at Bottom Level    
									cX1 = cX2 = grid->i2c[c_id].first;
									cY1 = cY2 = grid->i2c[c_id].second;
									e_ll[0] = cX1*(grid->cSide);        		
	    							e_ll[1] = cY1*(grid->cSide);        
									e_ur[0] = (cX2 + 1)*(grid->cSide);  
									e_ur[1] = (cY2 + 1)*(grid->cSide);
									priorityEntry.bottomLevel = 1;
								}								


								// Compute the zone
								if (priorityEntry.bottomLevel == 1) {  // ENTITY Is or Adjusted to Bottom Level
									e_center[0] = cX1 + 0.5;  // cell center (x) in terms of cell side			
									e_center[1] = cY1 + 0.5;  // cell center (y) in terms of cell side
									priorityEntry.zone = max(fabs(e_center[0] - q_center[0]), fabs(e_center[1] - q_center[1]));
									//cout << "zone: " << priorityEntry.zone << endl;	
								}	

								// COMPUTE THE 'MINDIST' OF RETRIEVED ENTITY < START >
								if (priorityEntry.zone == 0) {  // 'ZERO ZONE' case
									priorityEntry.minDist = 0.0;	
								}
								else {  // 'NON-ZERO ZONE' case
									// Check all the possible cases	
									if (e_ur[0] <= q_ll[0] && e_ur[1] <= q_ll[1]) {  // LOW-LEFT case
										// select the 'up-right' point (e_ur[0], e_ur[1])  
										p2[0] = (e_ur[0] * (grid->maxX - grid->minX)) + grid->minX;
										p2[1] = (e_ur[1] * (grid->maxY - grid->minY)) + grid->minY; 
										priorityEntry.minDist = distance2(p1, p2); 		
									}
									else if (e_ll[0] >= q_ur[0] && e_ur[1] <= q_ll[1]) {  // LOW-RIGHT case
										// select the 'up-left' point (e_ll[0], e_ur[1])
										p2[0] = (e_ll[0] * (grid->maxX - grid->minX)) + grid->minX;
										p2[1] = (e_ur[1] * (grid->maxY - grid->minY)) + grid->minY;
										priorityEntry.minDist = distance2(p1, p2);
									}  	
									else if (e_ur[0] <= q_ll[0] && e_ll[1] >= q_ur[1]) {  // UP-LEFT case
										// select the 'low-right' point (e_ur[0], e_ll[1])
										p2[0] = (e_ur[0] * (grid->maxX - grid->minX)) + grid->minX;
										p2[1] = (e_ll[1] * (grid->maxY - grid->minY)) + grid->minY;	
										priorityEntry.minDist = distance2(p1, p2);
									}
									else if (e_ll[0] >= q_ur[0] && e_ll[1] >= q_ur[1]) {  // UP-RIGHT case
										// select the 'low-left' point (e_ll[0], e_ll[1])
										p2[0] = (e_ll[0] * (grid->maxX - grid->minX)) + grid->minX;
										p2[1] = (e_ll[1] * (grid->maxY - grid->minY)) + grid->minY;	
										priorityEntry.minDist = distance2(p1, p2);
									}
									else if (e_ll[0] == q_ll[0] && e_ll[1] > q_ll[1]) {  // DIRECT-UP case
										p2[0] = nx1;
										p2[1] = q_ur[1] + ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'y' coord 
										p2[1] = (p2[1] * (grid->maxY - grid->minY)) + grid->minY;  // normalized 'y' coord
										priorityEntry.minDist = distance2(p1, p2);	
									}
									else if (e_ll[0] > q_ll[0] && e_ll[1] == q_ll[1]) {  // DIRECT-RIGHT case 
										p2[0] = q_ur[0] + ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'x' coord
										p2[0] = (p2[0] * (grid->maxX - grid->minX)) + grid->minX;  // normalized 'x' coord
										p2[1] = ny1;
										priorityEntry.minDist = distance2(p1, p2);
									}	
									else if (e_ll[0] == q_ll[0] && e_ll[1] < q_ll[1]) {  // DIRECT-DOWN case
										p2[0] = nx1;
										p2[1] = q_ll[1] - ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'y' coord
										p2[1] = (p2[1] * (grid->maxY - grid->minY)) + grid->minY;  // normalized 'y' coord	
										priorityEntry.minDist = distance2(p1, p2);
									}		
									else if (e_ll[0] < q_ll[0] && e_ll[1] == q_ll[1]) {  // DIRECT-LEFT case
										p2[0] = q_ll[0] - ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'x' coord
										p2[0] = (p2[0] * (grid->maxX - grid->minX)) + grid->minX;  // normalized 'x' coord		
										p2[1] = ny1;  
										priorityEntry.minDist = distance2(p1, p2);
									}
									else {
										cout << "MINDIST ERROR (3)!!" << endl;
										exit(-1);
									}			
								}
								// COMPUTE THE 'MINDIST' OF RETRIEVED ENTITY < END >

								//cout << "minDist: " << priorityEntry.minDist << endl;
								priorityEntry.dir = 4;  // cell
								H.push(priorityEntry);				
						
								// Stop taking more entities from previous operators ('spatial entity' case)
								if (subValue >= maxLimit) { 
									prevMaxLimit = subValue;
				   				break;
								}
							}
							else {  // All the entities are taken OR a non-spatial one is found that is above 'maxLimit'
				   			break;
							}
						}		
						// GET & PROCESS ALL THE SPATIAL ENTITIES FROM PREVIOUS OPERATORS, BASED ON THE  
						// '(priorityEntry.zone, priorityEntry.dir)' 'LIMIT ENTITY CONDITIONS' < END >
						}
					}		
				}
			}
			// APPLY CPM < END >

		
			// COLLECT KNN RESULTS < START >
			unsigned temp = 0;
			while (!R.empty()) {		
				priorityEntry = R.top();
				R.pop();
				kNN_buffer.push_back(kNN_res());
				kNN_buffer[temp].subject = priorityEntry.subject;
 				kNN_buffer[temp].geometry = priorityEntry.geometry; 
				kNN_buffer[temp].distance = priorityEntry.realDist;
				for (unsigned no_args = 0; no_args < priorityEntry.restArgs.size(); no_args++)
					kNN_buffer[temp].restArgs.push_back(priorityEntry.restArgs[no_args]);
				temp++; 
			}

			// Sort the buffer in ascending distance order
			sort(kNN_buffer.begin(), kNN_buffer.end(), sortByDistance());
			// COLLECT KNN RESULTS < END >


			// PASS REGISTER VALUES TO NEXT OPERATORS < START >
			cnt = 0;	 pos = 0;
			for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
				if ((*iter)->select_exist) {
					if ( (*iter)->subVar_exist ) 
						(*iter)->value = kNN_buffer[cnt].subject;
					else if ( (*iter)->geoVar_exist )
						(*iter)->value = kNN_buffer[cnt].geometry;
					else  
						(*iter)->value = kNN_buffer[cnt].restArgs[pos++]; 
				}
			}
			// PASS REGISTER VALUES TO NEXT OPERATORS < END >

  			observedOutputCardinality++;
  			cnt++;

  			return 1; 

		#else

			//cout << "first...kNN ENCODING UNSORTED Selection Evaluation" << endl;

			observedOutputCardinality = 0;			
			verif = filtered = 0;			

			// Get the first tuple
			unsigned count = entry->first();
			if (!count) return 0;
			
			// Get the 's' value of triple 's <hasGeometry> geo'
			unsigned subValue;
			for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
				if ( (*iter)->subVar_exist ) { 
					subValue = (*iter)->value;  // get and store the subject id 
					break;		
				}
			}


			// DEFINE 'QUERY' RELATED DATA < START >
			type1 = geo->type;
			//cout << "type1: " << type1 << endl;

			if (type1 == 1) {  // POINT coordinates in the Unit Square	   
   			rx1 = rx2 = ((Geo::Point*) geo)->getX(); 
   			ry1 = ry2 = ((Geo::Point*) geo)->getY(); 	
				// Coords Conversion from Unit Square to Normalized
				nx1 = nx2 = (rx1 * (grid->maxX - grid->minX)) + grid->minX; 	
				ny1 = ny2 = (ry1 * (grid->maxY - grid->minY)) + grid->minY;
			}
			else if (type1 == 2) {  // RECTANGLE coordinates in the Unit Square  
   			rx1 = ((Geo::Rectangle*) geo)->getP1()->getX(); 
   			ry1 = ((Geo::Rectangle*) geo)->getP1()->getY(); 
   			rx2 = ((Geo::Rectangle*) geo)->getP2()->getX(); 
   			ry2 = ((Geo::Rectangle*) geo)->getP2()->getY();		
				// Coords Conversion from Unit Square to Normalized
				nx1 = (rx1 * (grid->maxX - grid->minX)) + grid->minX;
				ny1 = (ry1 * (grid->maxY - grid->minY)) + grid->minY;
				nx2 = (rx2 * (grid->maxX - grid->minX)) + grid->minX;
				ny2 = (ry2 * (grid->maxY - grid->minY)) + grid->minY;
			}

			//cout << "(rx1, nx1): (" << rx1 << ", " << nx1 << ")" << endl;
			//cout << "(ry1, ny1): (" << ry1 << ", " << ny1 << ")" << endl;
			//cout << "(rx2, nx2): (" << rx2 << ", " << nx2 << ")" << endl;
			//cout << "(ry2, ny2): (" << ry2 << ", " << ny2 << ")" << endl;

			// get the cell id prefix "c_id" and the cell level "level"
			c_id = getHilbertId(rx1, rx2, ry1, ry2, grid->cSide, grid->c2i, grid->CELLS_PER_DIMENSION, grid->b1, grid->b2, level);
			// get the cell gpos "gpos" 
			gpos = getGridPos(c_id, level, grid->size, grid->MAX_LEVEL, grid->b1);

			if (!level) { 
				cout << "QUERY belongs to the Whole Grid" << endl;
				q_ll[0] = 0.0;												    
				q_ll[1] = 0.0;													  
				q_ur[0] = 1.0;  		
				q_ur[1] = 1.0;	 
				////////////////
				cout << "UNDEFINED YET ERROR (1)!!" << endl;
				exit(-1);
			}
			else if (level < grid->MAX_LEVEL) {
				cout << "QUERY is at Upper Level" << endl;
				gpos -= (((unsigned) 1) << grid->b1);
				cX1 = grid->offsets[gpos][0];
				cY1 = grid->offsets[gpos][1];
				cX2 = grid->offsets[gpos][2];
				cY2 = grid->offsets[gpos][3];	
				q_ll[0] = cX1*(grid->cSide);        		
		    	q_ll[1] = cY1*(grid->cSide);        
				q_ur[0] = (cX2 + 1)*(grid->cSide);  
				q_ur[1] = (cY2 + 1)*(grid->cSide);  
				////////////////
				cout << "UNDEFINED YET ERROR (2)!!" << endl;
				exit(-1);
			}
			else {
				//cout << "QUERY is at Bottom Level" << endl;
				cX1 = cX2 = grid->i2c[c_id].first;
				cY1 = cY2 = grid->i2c[c_id].second;
				q_ll[0] = cX1*(grid->cSide);        
		    	q_ll[1] = cY1*(grid->cSide);        
				q_ur[0] = (cX2 + 1)*(grid->cSide);  
				q_ur[1] = (cY2 + 1)*(grid->cSide);  

				q_center[0] = floor(q_ll[0] / grid->cSide) + 0.5;
				q_center[1] = floor(q_ll[1] / grid->cSide) + 0.5;
			}			

			//cout << "(Q[0], Q[1]): (" << q_center[0] << ", " << q_center[1] << ")" << endl << endl;
			// DEFINE 'QUERY' RELATED DATA < END >


			// INIT DIST < START >
			double p1[2];
			p1[0] = nx1;  
			p1[1] = ny1;  // QUERY coords

			double p2[2];
			p2[0] = nx1;
			p2[1] = (q_ur[1] * (grid->maxY - grid->minY)) + grid->minY;
			initDist[0] = distance2(p1, p2);  // U
			//cout << "initDist[0]: " << initDist[0] << endl;	
			p2[0] = (q_ur[0] * (grid->maxX - grid->minX)) + grid->minX;
			p2[1] = ny1;
			initDist[1] = distance2(p1, p2);  // R
			//cout << "initDist[1]: " << initDist[1] << endl;	  
			p2[0] = nx1;	
			p2[1] = (q_ll[1] * (grid->maxY - grid->minY)) + grid->minY;
			initDist[2] = distance2(p1, p2);  // D
			//cout << "initDist[2]: " << initDist[2] << endl;
			p2[0] = (q_ll[0] * (grid->maxX - grid->minX)) + grid->minX;		
			p2[1] = ny1;  
			initDist[3] = distance2(p1, p2);  // L
			//cout << "initDist[3]: " << initDist[3] << endl;			
			// INIT DIST < END >


			// GET & PROCESS THE FIRST SPATIAL ENTITY FROM PREVIOUS OPERATORS < START >			
			// Overlook the non-spatial entities
			while ( !getCellId(subValue, c_id, gpos, level) ) {  // get the cell "c_id", "gpos" and "level"
				entry->next();  // get the next tuple 
				// Get the 's' value of triple 's <hasGeometry> geo'
				for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
					if ( (*iter)->subVar_exist ) { 
						subValue = (*iter)->value;  // get and store the subject id
						break;
					} 		
				}
			}
	
			// Get and store the entity id
			priorityEntry.subject = subValue; 

			// Get and store the 'geo' id and all the 'rest variable' ids 
			for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
				if ( (*iter)->subVar_exist )
					continue; 		
				else if ( (*iter)->geoVar_exist )
					priorityEntry.geometry = (*iter)->value;  // get the ?geo id
				else  
					priorityEntry.restArgs.push_back((*iter)->value);
			}
	

			if (!level) {  // ENTITY belongs to the Whole Grid
				priorityEntry.zone = 0;
				priorityEntry.bottomLevel = 0;	
			}
			else if (level < grid->MAX_LEVEL) {  // ENTITY is at Upper Level
				// Take the REAL 'lower left' and 'upper right' points of the 'upper cell'
				gpos -= (((unsigned) 1) << grid->b1);
				cX1 = grid->offsets[gpos][0];
				cY1 = grid->offsets[gpos][1];
				cX2 = grid->offsets[gpos][2];
				cY2 = grid->offsets[gpos][3];	
				e_ll[0] = cX1*(grid->cSide);        		
		   	e_ll[1] = cY1*(grid->cSide);        
				e_ur[0] = (cX2 + 1)*(grid->cSide);  
				e_ur[1] = (cY2 + 1)*(grid->cSide);
  		
				// FIND THE NEAREST TO THE 'QUERY' BOTTOM CELL PREFIXES OF THE UPPER CELL < START >		
				if (q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // INTERSECTION case
					// select the 'query' cell
					cX1 = q_ll[0] / grid->cSide;  cY1 = q_ll[1] / grid->cSide;
				}
				else if (e_ur[0] <= q_ll[0] && e_ll[1] >= q_ur[1]) {  // UP-LEFT case
					// select the 'low-right' cell
					cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = e_ll[1] / grid->cSide;
				} 
				else if (e_ur[0] <= q_ll[0] && e_ur[1] <= q_ll[1]) {  // LOW-LEFT case
					// select the 'up-right' cell
					cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 
				}
				else if (e_ll[0] >= q_ur[0] && e_ll[1] >= q_ur[1]) {  // UP-RIGHT case
					// select the 'low-left' cell
					cX1 = e_ll[0] / grid->cSide;  cY1 = e_ll[1] / grid->cSide;
				}
				else if (e_ll[0] >= q_ur[0] && e_ur[1] <= q_ll[1]) {  // LOW-RIGHT case  
					// select the 'up-left' cell
					cX1 = e_ll[0] / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 			  	
				}
				else if (e_ur[0] <= q_ll[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // LEFT case
					// select the 'direct-right' cell
					cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = q_ll[1] / grid->cSide;      
				}
				else if (e_ll[0] >= q_ur[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // RIGHT case
					// select the 'direct-left' cell
					cX1 = e_ll[0] / grid->cSide;  cY1 = q_ll[1] / grid->cSide; 				
				}
				else if (e_ll[1] >= q_ur[1] && q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0]) {  // UP case
					// select the 'direct-low' cell
					cX1 = q_ll[0] / grid->cSide;  cY1 = e_ll[1] / grid->cSide;  
				}
				else if (e_ur[1] <= q_ll[1] && q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0]) {  // LOW case
					// select the 'direct-up' cell	
					cX1 = q_ll[0] / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 
				}
				else { 
					cout << "NEAREST BOTTOM CELL ERROR (1)!!" << endl;
					exit(-1);
				}	
				// FIND THE NEAREST TO THE 'QUERY' BOTTOM CELL PREFIXES OF THE UPPER CELL < END >
 
				// Take the ADJUSTED 'lower left' and 'upper right' points  
				// of the nearest to query 'bottom cell' of the 'upper cell'	
				e_ll[0] = cX1 * grid->cSide;
				e_ll[1] = cY1 * grid->cSide;
				e_ur[0] = (cX1 + 1) * grid->cSide;
				e_ur[1] = (cY1 + 1) * grid->cSide;
				// Consider this 'bottom cell' for the kNN evaluation 
				priorityEntry.bottomLevel = 1;
			}
			else {  // ENTITY is at Bottom Level
				cX1 = cX2 = grid->i2c[c_id].first;
				cY1 = cY2 = grid->i2c[c_id].second;
				e_ll[0] = cX1*(grid->cSide);        		
	    		e_ll[1] = cY1*(grid->cSide);        
				e_ur[0] = (cX2 + 1)*(grid->cSide);  
				e_ur[1] = (cY2 + 1)*(grid->cSide);  
				priorityEntry.bottomLevel = 1;
			}				


			// Compute the zone
			if (priorityEntry.bottomLevel == 1) {  // ENTITY Is or Adjusted to Bottom Level
				e_center[0] = cX1 + 0.5;  // cell center (x) in terms of cell side			
				e_center[1] = cY1 + 0.5;  // cell center (y) in terms of cell side
				priorityEntry.zone = max(fabs(e_center[0] - q_center[0]), fabs(e_center[1] - q_center[1]));
				//cout << "zone: " << priorityEntry.zone << endl;		
			} 

			// COMPUTE THE 'MINDIST' OF RETRIEVED ENTITY < START >
			if (priorityEntry.zone == 0) {  // 'ZERO ZONE' case			
				priorityEntry.minDist = 0.0;
			}
			else {  // 'NON-ZERO ZONE' case
				// Check all the possible cases	
				if (e_ur[0] <= q_ll[0] && e_ur[1] <= q_ll[1]) {  // LOW-LEFT case
					// select the 'up-right' point (e_ur[0], e_ur[1])  
					p2[0] = (e_ur[0] * (grid->maxX - grid->minX)) + grid->minX; 	
					p2[1] = (e_ur[1] * (grid->maxY - grid->minY)) + grid->minY;
					priorityEntry.minDist = distance2(p1, p2);	
				}
				else if (e_ll[0] >= q_ur[0] && e_ur[1] <= q_ll[1]) {  // LOW-RIGHT case
					// select the 'up-left' point (e_ll[0], e_ur[1])
					p2[0] = (e_ll[0] * (grid->maxX - grid->minX)) + grid->minX;	
					p2[1] = (e_ur[1] * (grid->maxY - grid->minY)) + grid->minY;
					priorityEntry.minDist = distance2(p1, p2);	
				}  	
				else if (e_ur[0] <= q_ll[0] && e_ll[1] >= q_ur[1]) {  // UP-LEFT case
					// select the 'low-right' point (e_ur[0], e_ll[1])
					p2[0] = (e_ur[0] * (grid->maxX - grid->minX)) + grid->minX; 		
					p2[1] = (e_ll[1] * (grid->maxY - grid->minY)) + grid->minY;
					priorityEntry.minDist = distance2(p1, p2); 	
				}
				else if (e_ll[0] >= q_ur[0] && e_ll[1] >= q_ur[1]) {  // UP-RIGHT case
					// select the 'low-left' point (e_ll[0], e_ll[1])
					p2[0] = (e_ll[0] * (grid->maxX - grid->minX)) + grid->minX;
					p2[1] = (e_ll[1] * (grid->maxY - grid->minY)) + grid->minY;	
					priorityEntry.minDist = distance2(p1, p2);	
				}
				else if (e_ll[0] == q_ll[0] && e_ll[1] > q_ll[1]) {  // DIRECT-UP case
					p2[0] = nx1;
					p2[1] = q_ur[1] + ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'y' coord 
					p2[1] = (p2[1] * (grid->maxY - grid->minY)) + grid->minY;  // normalized 'y' coord
					priorityEntry.minDist = distance2(p1, p2);	
				}
				else if (e_ll[0] > q_ll[0] && e_ll[1] == q_ll[1]) {  // DIRECT-RIGHT case 
					p2[0] = q_ur[0] + ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'x' coord
					p2[0] = (p2[0] * (grid->maxX - grid->minX)) + grid->minX;  // normalized 'x' coord
					p2[1] = ny1;
					priorityEntry.minDist = distance2(p1, p2);	
				}		
				else if (e_ll[0] == q_ll[0] && e_ll[1] < q_ll[1]) {  // DIRECT-DOWN case
					p2[0] = nx1;
					p2[1] = q_ll[1] - ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'y' coord
					p2[1] = (p2[1] * (grid->maxY - grid->minY)) + grid->minY;  // normalized 'y' coord	
					priorityEntry.minDist = distance2(p1, p2);					
				}	
				else if (e_ll[0] < q_ll[0] && e_ll[1] == q_ll[1]) {  // DIRECT-LEFT case
					p2[0] = q_ll[0] - ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'x' coord
					p2[0] = (p2[0] * (grid->maxX - grid->minX)) + grid->minX;  // normalized 'x' coord		
					p2[1] = ny1;  
					priorityEntry.minDist = distance2(p1, p2);					
				}
				else {
					cout << "MINDIST ERROR (1)!!" << endl;
					exit(-1);
				}
			} 			
			// COMPUTE THE 'MINDIST' OF RETRIEVED ENTITY < END >		

			//cout << "minDist: " << priorityEntry.minDist << endl;
			H.push(priorityEntry);	
			// GET & PROCESS THE FIRST SPATIAL ENTITY FROM PREVIOUS OPERATORS < END >


			// GET & PROCESS THE NEXT SPATIAL ENTITIES FROM PREVIOUS OPERATORS < START >
			while ( (count = entry->next()) != 0 ) {
				// Get the 's' value of triple 's <hasGeometry> geo'
				for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
					if ( (*iter)->subVar_exist ) { 
						subValue = (*iter)->value;  // get and store the subject id		
						break;
					} 		
				}

				// Overlook the non-spatial entities
				while ( !getCellId(subValue, c_id, gpos, level) ) {  // get the cell "c_id", "gpos" and "level"
					count = entry->next();  // get the next tuple
					if (!count) break;
					// Get the 's' value of triple 's <hasGeometry> geo'
					for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
						if ( (*iter)->subVar_exist ) { 
							subValue = (*iter)->value;  // get and store the subject id
							break;
						} 		
					}
				}

				if (count != 0) {  // Edit all the spatial entities

					// Clear previous 'argument values'
					priorityEntry.restArgs.clear();

					// Get and store the entity id
					priorityEntry.subject = subValue; 
	
					// Get and store the 'geo' id and all the 'rest variable' ids 
					for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
						if ( (*iter)->subVar_exist )
							continue; 		
						else if ( (*iter)->geoVar_exist )
							priorityEntry.geometry = (*iter)->value;  // get the ?geo id
						else 
							priorityEntry.restArgs.push_back((*iter)->value);
					}

					
					if (!level) {  // ENTITY belongs to the Whole Grid  
						priorityEntry.zone = 0; 	
						priorityEntry.bottomLevel = 0; 
					}
					else if (level < grid->MAX_LEVEL) {  // ENTITY is at Upper Level
						// Take the REAL 'lower left' and 'upper right' points of the 'upper cell'
						gpos -= (((unsigned) 1) << grid->b1);
						cX1 = grid->offsets[gpos][0];
						cY1 = grid->offsets[gpos][1];
						cX2 = grid->offsets[gpos][2];
						cY2 = grid->offsets[gpos][3];	
						e_ll[0] = cX1*(grid->cSide);       		
		   			e_ll[1] = cY1*(grid->cSide);        
						e_ur[0] = (cX2 + 1)*(grid->cSide); 
						e_ur[1] = (cY2 + 1)*(grid->cSide); 

						// FIND THE NEAREST TO THE 'QUERY' BOTTOM CELL PREFIXES OF THE UPPER CELL < START >		
						if (q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // INTERSECTION case
							// select the 'query' cell
							cX1 = q_ll[0] / grid->cSide;  cY1 = q_ll[1] / grid->cSide;
						}
						else if (e_ur[0] <= q_ll[0] && e_ll[1] >= q_ur[1]) {  // UP-LEFT case
							// select the 'low-right' cell
							cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = e_ll[1] / grid->cSide;
						} 
						else if (e_ur[0] <= q_ll[0] && e_ur[1] <= q_ll[1]) {  // LOW-LEFT case
							// select the 'up-right' cell
							cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 
						}
						else if (e_ll[0] >= q_ur[0] && e_ll[1] >= q_ur[1]) {  // UP-RIGHT case
							// select the 'low-left' cell
							cX1 = e_ll[0] / grid->cSide;  cY1 = e_ll[1] / grid->cSide;
						}
						else if (e_ll[0] >= q_ur[0] && e_ur[1] <= q_ll[1]) {  // LOW-RIGHT case  
							// select the 'up-left' cell
							cX1 = e_ll[0] / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 			  	
						}
						else if (e_ur[0] <= q_ll[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // LEFT case
							// select the 'direct-right' cell
							cX1 = (e_ur[0] - grid->cSide) / grid->cSide;  cY1 = q_ll[1] / grid->cSide;      
						}
						else if (e_ll[0] >= q_ur[0] && q_ll[1] >= e_ll[1] && q_ur[1] <= e_ur[1]) {  // RIGHT case
							// select the 'direct-left' cell
							cX1 = e_ll[0] / grid->cSide;  cY1 = q_ll[1] / grid->cSide; 				
						}
						else if (e_ll[1] >= q_ur[1] && q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0]) {  // UP case
							// select the 'direct-low' cell
							cX1 = q_ll[0] / grid->cSide;  cY1 = e_ll[1] / grid->cSide;  
						}
						else if (e_ur[1] <= q_ll[1] && q_ll[0] >= e_ll[0] && q_ur[0] <= e_ur[0]) {  // LOW case
							// select the 'direct-up' cell	
							cX1 = q_ll[0] / grid->cSide;  cY1 = (e_ur[1] - grid->cSide) / grid->cSide; 
						}
						else { 
							cout << "NEAREST BOTTOM CELL ERROR (2)!!" << endl;
							exit(-1);
						}	
						// FIND THE NEAREST TO THE 'QUERY' BOTTOM CELL PREFIXES OF THE UPPER CELL < END >

						// Take the ADJUSTED 'lower left' and 'upper right' points  
						// of the nearest to query 'bottom cell' of the 'upper cell'	
						e_ll[0] = cX1 * grid->cSide;
						e_ll[1] = cY1 * grid->cSide;
						e_ur[0] = (cX1 + 1) * grid->cSide;
						e_ur[1] = (cY1 + 1) * grid->cSide;
						// Consider this 'bottom cell' for the kNN evaluation 
						priorityEntry.bottomLevel = 1;
					}
					else {  // ENTITY is at Bottom Level    
						cX1 = cX2 = grid->i2c[c_id].first;
						cY1 = cY2 = grid->i2c[c_id].second;
						e_ll[0] = cX1*(grid->cSide);        		
	    				e_ll[1] = cY1*(grid->cSide);        
						e_ur[0] = (cX2 + 1)*(grid->cSide);  
						e_ur[1] = (cY2 + 1)*(grid->cSide);  
						priorityEntry.bottomLevel = 1;
					}								


					// Compute the zone
					if (priorityEntry.bottomLevel == 1) {  // ENTITY Is or Adjusted to Bottom Level
						e_center[0] = cX1 + 0.5;  // cell center (x) in terms of cell side			
						e_center[1] = cY1 + 0.5;  // cell center (y) in terms of cell side
						priorityEntry.zone = max(fabs(e_center[0] - q_center[0]), fabs(e_center[1] - q_center[1]));
						//cout << "zone: " << priorityEntry.zone << endl;	
					}	

					// COMPUTE THE 'MINDIST' OF RETRIEVED ENTITY < START >
					if (priorityEntry.zone == 0) {  // 'ZERO ZONE' case
						priorityEntry.minDist = 0.0;	
					}
					else {  // 'NON-ZERO ZONE' case
						// Check all the possible cases	
						if (e_ur[0] <= q_ll[0] && e_ur[1] <= q_ll[1]) {  // LOW-LEFT case
							// select the 'up-right' point (e_ur[0], e_ur[1])  
							p2[0] = (e_ur[0] * (grid->maxX - grid->minX)) + grid->minX;
							p2[1] = (e_ur[1] * (grid->maxY - grid->minY)) + grid->minY; 
							priorityEntry.minDist = distance2(p1, p2); 
						}
						else if (e_ll[0] >= q_ur[0] && e_ur[1] <= q_ll[1]) {  // LOW-RIGHT case
							// select the 'up-left' point (e_ll[0], e_ur[1])
							p2[0] = (e_ll[0] * (grid->maxX - grid->minX)) + grid->minX;
							p2[1] = (e_ur[1] * (grid->maxY - grid->minY)) + grid->minY;
							priorityEntry.minDist = distance2(p1, p2);	
						}  	
						else if (e_ur[0] <= q_ll[0] && e_ll[1] >= q_ur[1]) {  // UP-LEFT case
							// select the 'low-right' point (e_ur[0], e_ll[1])
							p2[0] = (e_ur[0] * (grid->maxX - grid->minX)) + grid->minX;
							p2[1] = (e_ll[1] * (grid->maxY - grid->minY)) + grid->minY;	
							priorityEntry.minDist = distance2(p1, p2); 	
						}
						else if (e_ll[0] >= q_ur[0] && e_ll[1] >= q_ur[1]) {  // UP-RIGHT case
							// select the 'low-left' point (e_ll[0], e_ll[1])
							p2[0] = (e_ll[0] * (grid->maxX - grid->minX)) + grid->minX;
							p2[1] = (e_ll[1] * (grid->maxY - grid->minY)) + grid->minY;	
							priorityEntry.minDist = distance2(p1, p2);	
						}
						else if (e_ll[0] == q_ll[0] && e_ll[1] > q_ll[1]) {  // DIRECT-UP case
							p2[0] = nx1;
							p2[1] = q_ur[1] + ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'y' coord 
							p2[1] = (p2[1] * (grid->maxY - grid->minY)) + grid->minY;  // normalized 'y' coord
							priorityEntry.minDist = distance2(p1, p2);	
						}
						else if (e_ll[0] > q_ll[0] && e_ll[1] == q_ll[1]) {  // DIRECT-RIGHT case 
							p2[0] = q_ur[0] + ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'x' coord
							p2[0] = (p2[0] * (grid->maxX - grid->minX)) + grid->minX;  // normalized 'x' coord
							p2[1] = ny1;
							priorityEntry.minDist = distance2(p1, p2);							
						}	
						else if (e_ll[0] == q_ll[0] && e_ll[1] < q_ll[1]) {  // DIRECT-DOWN case
							p2[0] = nx1;
							p2[1] = q_ll[1] - ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'y' coord
							p2[1] = (p2[1] * (grid->maxY - grid->minY)) + grid->minY;  // normalized 'y' coord	
							priorityEntry.minDist = distance2(p1, p2);							
						}	
						else if (e_ll[0] < q_ll[0] && e_ll[1] == q_ll[1]) {  // DIRECT-LEFT case
							p2[0] = q_ll[0] - ((priorityEntry.zone - 1) * grid->cSide);  // unit square 'x' coord
							p2[0] = (p2[0] * (grid->maxX - grid->minX)) + grid->minX;  // normalized 'x' coord		
							p2[1] = ny1;  
							priorityEntry.minDist = distance2(p1, p2);							
						}
						else {
							cout << "MINDIST ERROR (2)!!" << endl;
							exit(-1);
						}			
					}
					// COMPUTE THE 'MINDIST' OF RETRIEVED ENTITY < END >

					//cout << "minDist: " << priorityEntry.minDist << endl;
					H.push(priorityEntry);	
				}
				else {  // All the entities are taken
				   break;
				}
			}			
			// GET & PROCESS THE NEXT SPATIAL ENTITIES FROM PREVIOUS OPERATORS < END >


			// QUERY COORDS INITIALIZATION FOR DISTANCE FUNCTIONS < START >
			if (type1 == 1) {  // query is POINT
				size1 = 2;
				points1 = new double[size1];
				points1[0] = nx1;  points1[1] = ny1;
				type1 = 0;  // adjusted for DISTANCE functions
			}	
			else if (type1 == 2) {  // query is RECTANGLE
				size1 = 4;
				points1 = new double[size1];
				points1[0] = nx1;  points1[1] = ny1;
				points1[2] = nx2;  points1[3] = ny2;				
			}
			// QUERY COORDS INITIALIZATION FOR DISTANCE FUNCTIONS < END >


			// APPLY CPM < START >
			Result geoRes;
			bool f;
			// Only 'CELL ENTITY' cases exist in 'H' (Naive CPM Approach)
			while (!H.empty()) { 
				priorityEntry = H.top();
				if (priorityEntry.minDist >= lastDist)
					break; 
				// Remove the retrieved Entity from 'H'
				H.pop();  

            geoRes.flags = 0;  // reset the flags
				geoRes.flags = Result::idAvailable;
				geoRes.id = priorityEntry.geometry;  // get the ?geo id
				geoRes.ensureString(this);  // get the ?geo string value ('retrieve geometry')

				type2 = getLongLat(geoRes.value, points2, size2);  // set 'points2' & 'size2' values (NORMALIZED COORDS)

				// COMPUTE REAL DISTANCE < START >
				f = false;
		
				if (!type1 && !type2) {  // both are points
					//cout << "Point vs Point" << endl;
					priorityEntry.realDist = distance2(points1, points2);
				}
				else if ((f=(!type1 && type2==1)) || (!type2 && type1==1)) {  // point line
					if (f) {  // point line
						//cout << "Point vs Line" << endl;
						priorityEntry.realDist = pointLineStringDist2(points1, points2, size2);
					}
					else {  // line point
						//cout << "Line vs Point" << endl;
						priorityEntry.realDist = pointLineStringDist2(points2, points1, size1);
					}
				}
   			else if ((f=(!type1 && type2==2)) || (!type2 && type1==2)) {  // point polygon
					if (f) {  // point polygon
						//cout << "Point vs Polygon" << endl;
						priorityEntry.realDist = pointPolygonDist2(points1, points2, size2);
					}
					else {  // polygon point
						//cout << "Polygon vs Point" << endl;
						priorityEntry.realDist = pointPolygonDist2(points2, points1, size1);
					}	
   			}
				else if ((f=(type1==1 && type2==2)) || (type2==1 && type1==2)) {  // line polygon
					if (f) {  // line polygon
						//cout << "Line vs Polygon" << endl;
						priorityEntry.realDist = lineStringPolygonDist2(points1, size1, points2, size2);
					}
					else {  // polygon line
						//cout << "Polygon vs Line" << endl;
						priorityEntry.realDist = lineStringPolygonDist2(points2, size2, points1, size1);
					}
				}	
				else if (type1==1 && type2==1) {  // line line
					//cout << "Line vs Line" << endl;	
					priorityEntry.realDist = lineStringLineStringDist2(points1, size1, points2, size2);
				}
				else if (type1==2 && type2==2) {  // polygon polygon
					//cout << "Polygon vs Polygon" << endl;
					priorityEntry.realDist = polygonPolygonDist2(points1, points2, size1, size2);
				}
				else if ((f=(!type1 && type2==3)) || (!type2 && type1==3)) {  // point multipoint
					if (f) {  // point multipoint	
						//cout << "Point vs Multipoint" << endl;
						priorityEntry.realDist = pointMultipointDist2(points1, points2, size2);
					}
					else {  // multipoint point
						//cout << "Multipoint vs Point" << endl;
						priorityEntry.realDist = pointMultipointDist2(points2, points1, size1);
					}
				}
				else if (type1==3 && type2==3) {  // multipoint multipoint
					//cout << "Multipoint vs Multipoint" << endl;
					priorityEntry.realDist = multipointMultipointDist2(points1, points2, size1, size2);	
				}
				else if ((f=(type1==1 && type2==3)) || (type2==1 && type1==3)) {  // line multipoint
					if (f) {  // line multipoint
						//cout << "Line vs Multipoint" << endl;
						priorityEntry.realDist = lineStringMultipointDist2(points1, size1, points2, size2);
					}
					else {  // multipoint line
						//cout << "Multipoint vs Line" << endl;
						priorityEntry.realDist = lineStringMultipointDist2(points2, size2, points1, size1);
					}
				}
				else if ((f=(type1==2 && type2==3)) || (type2==2 && type1==3)) {  // polygon multipoint
					if (f) {  // polygon multipoint
						//cout << "Polygon vs Multipoint" << endl;
						priorityEntry.realDist = polygonMultipointDist2(points1, points2, size1, size2);
					}
					else {  // multipoint polygon
						//cout << "Multipoint vs Polygon" << endl;
						priorityEntry.realDist = polygonMultipointDist2(points2, points1, size2, size1);	
					}
				}
				else {
					cout << "kNN Selection ERROR: Could not recognize pair of geometries: " << type1 << " " << type2 << endl;
					exit(-1);
				}
				// COMPUTE REAL DISTANCE < END >


				cell_cnt++;
				if (cell_cnt > k_nearest_neighbor) {
					R.push(priorityEntry);				
					R.pop();
				}
				else {  
					R.push(priorityEntry);
				}

				if (cell_cnt >= k_nearest_neighbor) { 
					priorityEntry = R.top();
					lastDist = priorityEntry.realDist;
				}
			}
			// APPLY CPM < END >


			// COLLECT KNN RESULTS < START >
			unsigned temp = 0;
			while (!R.empty()) {		
				priorityEntry = R.top();
				R.pop();
				kNN_buffer.push_back(kNN_res());
				kNN_buffer[temp].subject = priorityEntry.subject;
 				kNN_buffer[temp].geometry = priorityEntry.geometry; 
				kNN_buffer[temp].distance = priorityEntry.realDist;
				for (unsigned no_args = 0; no_args < priorityEntry.restArgs.size(); no_args++)
					kNN_buffer[temp].restArgs.push_back(priorityEntry.restArgs[no_args]);
				temp++; 
			}

			// Sort the buffer in ascending distance order
			sort(kNN_buffer.begin(), kNN_buffer.end(), sortByDistance());
			// COLLECT KNN RESULTS < END >

		
			// PASS REGISTER VALUES TO NEXT OPERATORS < START >
			cnt = 0;  pos = 0;
			for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
				if ((*iter)->select_exist) {
					if ( (*iter)->subVar_exist ) 
						(*iter)->value = kNN_buffer[cnt].subject;
					else if ( (*iter)->geoVar_exist )
						(*iter)->value = kNN_buffer[cnt].geometry;
					else  
						(*iter)->value = kNN_buffer[cnt].restArgs[pos++]; 
				}
			}
			// PASS REGISTER VALUES TO NEXT OPERATORS < END >

  			observedOutputCardinality++;
  			cnt++;

  			return 1;

		#endif
	}

	else {

		#if RTREE_USED_IN_PLAN

			#if ONLY_POINT_GEOMS

				//cout << "first...kNN RTREE Selection Evaluation (POINT_GEOMS)" << endl;

				observedOutputCardinality = 0;
				verif = filtered = 0;

				// Get the first tuple
				unsigned count = entry->first();
				if (!count) return 0;


				kNN_buffer.push_back(kNN_res());  // allocate a new entry to kNN buffer

				for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
					if ( (*iter)->subVar_exist ) 
						kNN_buffer[cnt].subject = (*iter)->value;  // get and store the subject id 		
					else if ( (*iter)->geoVar_exist )
						continue;
					else 
						kNN_buffer[cnt].restArgs.push_back((*iter)->value);
				}
				cnt++;

      		// Get the next tuples
				while ((cnt != k_nearest_neighbor) && (entry->next() != 0)) {

					kNN_buffer.push_back(kNN_res());  // allocate a new entry to kNN buffer

					for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
						if ( (*iter)->subVar_exist ) 
							kNN_buffer[cnt].subject = (*iter)->value;  // get and store the subject id 		
						else if ( (*iter)->geoVar_exist )
							continue;	
						else 
							kNN_buffer[cnt].restArgs.push_back((*iter)->value);	
					}
					cnt++;
				}


				// PASS REGISTER VALUES TO NEXT OPERATORS < START >
				cnt = 0;  pos = 0;
				for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
					if ((*iter)->select_exist) {
						if ( (*iter)->subVar_exist ) 
							(*iter)->value = kNN_buffer[cnt].subject;
						else if ( (*iter)->geoVar_exist )
							continue;
						else  
							(*iter)->value = kNN_buffer[cnt].restArgs[pos++]; 
					}
				}
				// PASS REGISTER VALUES TO NEXT OPERATORS < END >

  				observedOutputCardinality++;
  				cnt++;

  				return 1;

			#else

				//cout << "first...kNN RTREE Selection Evaluation (ALL_GEOMS)" << endl;

				observedOutputCardinality = 0;
				verif = filtered = 0;

				// Get the first tuple
				unsigned count = entry->first();
				if (!count) return 0;


				// DEFINE 'QUERY' RELATED DATA < START >
				type1 = geo->type;	
				//cout << "type1: " << type1 << endl;	
				if (type1 == 1) {  // query is POINT
					size1 = 2;
					points1 = new double[size1];
					// Take Normalized Coords
					points1[0] = ((Geo::Point*) geo)->getX();
					points1[1] = ((Geo::Point*) geo)->getY();  	
					//cout << "nX_1: " << points1[0] << endl;  cout << "nY_1: " << points1[1] << endl;
					type1 = 0;  // adjusted for DISTANCE functions	
				}
				else if (type1 == 2) {  // query is RECTANGLE  
   				size1 = 4;			
					points1 = new double[size1];
					// Take Normalized Coords	
					points1[0] = ((Geo::Rectangle*) geo)->getP1()->getX(); 
   				points1[1] = ((Geo::Rectangle*) geo)->getP1()->getY(); 
   				points1[2] = ((Geo::Rectangle*) geo)->getP2()->getX(); 
   				points1[3] = ((Geo::Rectangle*) geo)->getP2()->getY();
					//cout << "nX_1: " << points1[0] << endl;  cout << "nY_1: " << points1[1] << endl;
					//cout << "nX_2: " << points1[2] << endl;  cout << "nY_2: " << points1[3] << endl;
				}
				// DEFINE 'QUERY' RELATED DATA < END >


				// VARIABLE DECLARATION < START >
				#if LGD
					string hasGeom_string = "hasGeometry";
				#elif YAGO
					string hasGeom_string = "http://yago-knowledge.org/resource/hasGeometry";
				#endif

				Type::ID hasGeom_type;
				unsigned hasGeom_subType;
				unsigned hasGeom_id; 
				bool exist;

				exist = runtime.getDatabase().getDictionary().lookup(hasGeom_string, hasGeom_type, hasGeom_subType, hasGeom_id);
				if (exist) {
					//cout << "HASGEO_ID: " << hasGeom_id << endl;
				}
				else { 
					cout << "HASGEO_ID NOT FOUND ERROR!" << endl;
					exit(-1);
				}
	
				unsigned filter1;  // the ?s variable value
				unsigned filter2 = hasGeom_id;  // the <hasGeometry> constant value 
				Result geoRes;
				bool f;
				double minX, maxX, minY, maxY;
				double MBR_points[4];
				double MBR_dist;
				// VARIABLE DECLARATION < END >


				// Process the first tuple
				for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
					if ( (*iter)->subVar_exist ) 
						prRtreeEntry.subject = (*iter)->value;  // get the subject id 		
					else if ( (*iter)->geoVar_exist )
						continue;
					else 
						prRtreeEntry.restArgs.push_back((*iter)->value);
				} 
				filter1 = prRtreeEntry.subject;
				prRtreeEntry.geometry = IS->getGeoID(filter1, filter2);  // get the ?geo id

				geoRes.flags = 0;  // reset the flags
				geoRes.flags = Result::idAvailable;
  				geoRes.id = prRtreeEntry.geometry;
				geoRes.ensureString(this);  // get the ?geo string value ('retrieve geometry')  

				#if YAGO
					spatialIDs_TO_geomIDs.insert(make_pair(prRtreeEntry.subject, prRtreeEntry.geometry));
				#endif	

				type2 = getLongLat(geoRes.value, points2, size2);  // set 'points2' & 'size2' values (NORMALIZED COORDS)

				// COMPUTE REAL DISTANCE < START >
				f = false;

				if (!type1 && !type2) {  // both are points
					//cout << "Point vs Point" << endl;
					prRtreeEntry.distance = distance2(points1, points2);
				}
				else if ((f=(!type1 && type2==1)) || (!type2 && type1==1)) {  // point line
					if (f) {  // point line
						//cout << "Point vs Line" << endl;
						prRtreeEntry.distance = pointLineStringDist2(points1, points2, size2);
					}
					else {  // line point
						//cout << "Line vs Point" << endl;
						prRtreeEntry.distance = pointLineStringDist2(points2, points1, size1);
					}
				}
   			else if ((f=(!type1 && type2==2)) || (!type2 && type1==2)) {  // point polygon
					if (f) {  // point polygon
						//cout << "Point vs Polygon" << endl;
						prRtreeEntry.distance = pointPolygonDist2(points1, points2, size2);
					}
					else {  // polygon point
						//cout << "Polygon vs Point" << endl;
						prRtreeEntry.distance = pointPolygonDist2(points2, points1, size1);
					}	
   			}
				else if ((f=(type1==1 && type2==2)) || (type2==1 && type1==2)) {  // line polygon
					if (f) {  // line polygon
						//cout << "Line vs Polygon" << endl;
						prRtreeEntry.distance = lineStringPolygonDist2(points1, size1, points2, size2);
					}
					else {  // polygon line
						//cout << "Polygon vs Line" << endl;
						prRtreeEntry.distance = lineStringPolygonDist2(points2, size2, points1, size1);
					}
				}
				else if (type1==1 && type2==1) {  // line line
					//cout << "Line vs Line" << endl;	
					prRtreeEntry.distance = lineStringLineStringDist2(points1, size1, points2, size2);
				}
				else if (type1==2 && type2==2) {  // polygon polygon
					//cout << "Polygon vs Polygon" << endl;
					prRtreeEntry.distance = polygonPolygonDist2(points1, points2, size1, size2);
				}
				else if ((f=(!type1 && type2==3)) || (!type2 && type1==3)) {  // point multipoint
					if (f) {  // point multipoint	
						//cout << "Point vs Multipoint" << endl;
						prRtreeEntry.distance = pointMultipointDist2(points1, points2, size2);
					}
					else {  // multipoint point
						//cout << "Multipoint vs Point" << endl;
						prRtreeEntry.distance = pointMultipointDist2(points2, points1, size1);
					}
				}
				else if (type1==3 && type2==3) {  // multipoint multipoint
					//cout << "Multipoint vs Multipoint" << endl;
					prRtreeEntry.distance = multipointMultipointDist2(points1, points2, size1, size2);	
				}
				else if ((f=(type1==1 && type2==3)) || (type2==1 && type1==3)) {  // line multipoint
					if (f) {  // line multipoint
						//cout << "Line vs Multipoint" << endl;
						prRtreeEntry.distance = lineStringMultipointDist2(points1, size1, points2, size2);
					}
					else {  // multipoint line
						//cout << "Multipoint vs Line" << endl;
						prRtreeEntry.distance = lineStringMultipointDist2(points2, size2, points1, size1);
					}
				}
				else if ((f=(type1==2 && type2==3)) || (type2==2 && type1==3)) {  // polygon multipoint
					if (f) {  // polygon multipoint
						//cout << "Polygon vs Multipoint" << endl;
						prRtreeEntry.distance = polygonMultipointDist2(points1, points2, size1, size2);
					}
					else {  // multipoint polygon
						//cout << "Multipoint vs Polygon" << endl;
						prRtreeEntry.distance = polygonMultipointDist2(points2, points1, size2, size1);	
					}
				}
				else {
					cout << "kNN Selection ERROR: Could not recognize pair of geometries: " << type1 << " " << type2 << endl;
					exit(-1);
				}
				// COMPUTE REAL DISTANCE < END >

				RR.push(prRtreeEntry);
				cnt++;
		

      		// Get & Process the next tuples
				while (entry->next() != 0) {

					// Clear previous 'argument values'
					prRtreeEntry.restArgs.clear();

					for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
						if ( (*iter)->subVar_exist ) 
							prRtreeEntry.subject = (*iter)->value;  // get the subject id	 		
						else if ( (*iter)->geoVar_exist )
							continue;	
						else 
							prRtreeEntry.restArgs.push_back((*iter)->value);	
					}

					geoRes.flags = 0;  // reset the flags
					geoRes.flags = Result::idAvailable;

					#if LGD

						filter1 = prRtreeEntry.subject;
						prRtreeEntry.geometry = IS->getGeoID(filter1, filter2);  // get the ?geo id
  						geoRes.id = prRtreeEntry.geometry;

  					#elif YAGO	

						auto find_it = spatialIDs_TO_geomIDs.find(prRtreeEntry.subject);

						if (find_it == spatialIDs_TO_geomIDs.end()) {  // not found
							filter1 = prRtreeEntry.subject;
							prRtreeEntry.geometry = IS->getGeoID(filter1, filter2);  // get the ?geo id
							geoRes.id = prRtreeEntry.geometry;
							spatialIDs_TO_geomIDs.insert(make_pair(prRtreeEntry.subject, prRtreeEntry.geometry));	
						}	
						else {  // found   
							geoRes.id = find_it->second; 
						}

  					#endif	

					geoRes.ensureString(this);  // get the ?geo string value ('retrieve geometry')  

					type2 = getLongLat(geoRes.value, points2, size2);  // set 'points2' & 'size2' values (NORMALIZED COORDS)

					getRegion(points2, size2, minX, maxX, minY, maxY);  // get the MBR of the geometry	
					MBR_points[0] = minX;  MBR_points[1] = minY;  MBR_points[2] = maxX;  MBR_points[3] = maxY; 								
					MBR_dist = RtreeDistPointRegion(points1, MBR_points, 2);

					// COMPUTE REAL DISTANCE < START >
					f = false;

					if (!type1 && !type2) {  // both are points
						//cout << "Point vs Point" << endl;
						prRtreeEntry.distance = distance2(points1, points2);
					}
					else if ((f=(!type1 && type2==1)) || (!type2 && type1==1)) {  // point line
						if (f) {  // point line
							//cout << "Point vs Line" << endl;
							prRtreeEntry.distance = pointLineStringDist2(points1, points2, size2);
						}
						else {  // line point
							//cout << "Line vs Point" << endl;
							prRtreeEntry.distance = pointLineStringDist2(points2, points1, size1);
						}
					}
   				else if ((f=(!type1 && type2==2)) || (!type2 && type1==2)) {  // point polygon
						if (f) {  // point polygon
							//cout << "Point vs Polygon" << endl;
							prRtreeEntry.distance = pointPolygonDist2(points1, points2, size2);
						}
						else {  // polygon point
							//cout << "Polygon vs Point" << endl;
							prRtreeEntry.distance = pointPolygonDist2(points2, points1, size1);
						}	
   				}
					else if ((f=(type1==1 && type2==2)) || (type2==1 && type1==2)) {  // line polygon
						if (f) {  // line polygon
							//cout << "Line vs Polygon" << endl;
							prRtreeEntry.distance = lineStringPolygonDist2(points1, size1, points2, size2);
						}
						else {  // polygon line
							//cout << "Polygon vs Line" << endl;
							prRtreeEntry.distance = lineStringPolygonDist2(points2, size2, points1, size1);
						}
					}
					else if (type1==1 && type2==1) {  // line line
						//cout << "Line vs Line" << endl;	
						prRtreeEntry.distance = lineStringLineStringDist2(points1, size1, points2, size2);
					}
					else if (type1==2 && type2==2) {  // polygon polygon
						//cout << "Polygon vs Polygon" << endl;
						prRtreeEntry.distance = polygonPolygonDist2(points1, points2, size1, size2);
					}
					else if ((f=(!type1 && type2==3)) || (!type2 && type1==3)) {  // point multipoint
						if (f) {  // point multipoint	
							//cout << "Point vs Multipoint" << endl;
							prRtreeEntry.distance = pointMultipointDist2(points1, points2, size2);
						}
						else {  // multipoint point
							//cout << "Multipoint vs Point" << endl;
							prRtreeEntry.distance = pointMultipointDist2(points2, points1, size1);
						}
					}
					else if (type1==3 && type2==3) {  // multipoint multipoint
						//cout << "Multipoint vs Multipoint" << endl;
						prRtreeEntry.distance = multipointMultipointDist2(points1, points2, size1, size2);	
					}
					else if ((f=(type1==1 && type2==3)) || (type2==1 && type1==3)) {  // line multipoint
						if (f) {  // line multipoint
							//cout << "Line vs Multipoint" << endl;
							prRtreeEntry.distance = lineStringMultipointDist2(points1, size1, points2, size2);
						}
						else {  // multipoint line
							//cout << "Multipoint vs Line" << endl;
							prRtreeEntry.distance = lineStringMultipointDist2(points2, size2, points1, size1);
						}
					}
					else if ((f=(type1==2 && type2==3)) || (type2==2 && type1==3)) {  // polygon multipoint
						if (f) {  // polygon multipoint
							//cout << "Polygon vs Multipoint" << endl;
							prRtreeEntry.distance = polygonMultipointDist2(points1, points2, size1, size2);
						}
						else {  // multipoint polygon
							//cout << "Multipoint vs Polygon" << endl;
							prRtreeEntry.distance = polygonMultipointDist2(points2, points1, size2, size1);	
						}
					}
					else {
						cout << "kNN Selection ERROR: Could not recognize pair of geometries: " << type1 << " " << type2 << endl;
						exit(-1);
					}
					// COMPUTE REAL DISTANCE < END >
		
					if (cnt >= k_nearest_neighbor) {
						prRtreeEntryTop = RR.top();
						if (MBR_dist >= prRtreeEntryTop.distance)
							break;
						else {
							RR.push(prRtreeEntry);
							RR.pop();
						}	
					}
					else
						RR.push(prRtreeEntry);

					cnt++;
				}


				// COLLECT KNN RESULTS < START >
				unsigned temp = 0;
				while (!RR.empty()) {		
					prRtreeEntry = RR.top();
					RR.pop();
					kNN_buffer.push_back(kNN_res());
					kNN_buffer[temp].subject = prRtreeEntry.subject;
 					kNN_buffer[temp].geometry = prRtreeEntry.geometry; 
					kNN_buffer[temp].distance = prRtreeEntry.distance;
					for (unsigned no_args = 0; no_args < prRtreeEntry.restArgs.size(); no_args++)
						kNN_buffer[temp].restArgs.push_back(prRtreeEntry.restArgs[no_args]);
					temp++; 
				}

				// Sort the buffer in ascending distance order
				sort(kNN_buffer.begin(), kNN_buffer.end(), sortByDistance());
				// COLLECT KNN RESULTS < END >


				// PASS REGISTER VALUES TO NEXT OPERATORS < START >
				cnt = 0;  pos = 0;
				for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
					if ((*iter)->select_exist) {
						if ( (*iter)->subVar_exist ) 
							(*iter)->value = kNN_buffer[cnt].subject;
						else if ( (*iter)->geoVar_exist )
							(*iter)->value = kNN_buffer[cnt].geometry;
						else  
							(*iter)->value = kNN_buffer[cnt].restArgs[pos++]; 
					}
				}
				// PASS REGISTER VALUES TO NEXT OPERATORS < END >

  				observedOutputCardinality++;
  				cnt++;

  				return 1; 

			#endif

  		#else		

			cout.precision(10);
		
			//cout << "first...kNN BASELINE Selection Evaluation" << endl;
	
			observedOutputCardinality = 0;
			verif = filtered = 0;

			// Get the first tuple
			unsigned count = entry->first();
			if (!count) return 0;

			
			// DEFINE 'QUERY' RELATED DATA < START >
			type1 = geo->type;
			//cout << "type1: " << type1 << endl;
			if (type1 == 1) {  // query is POINT
				size1 = 2;
				points1 = new double[size1];
				// Take Normalized Coords
				points1[0] = ((Geo::Point*) geo)->getX();
				points1[1] = ((Geo::Point*) geo)->getY(); 	
				//cout << "nX_1: " << points1[0] << endl;  cout << "nY_1: " << points1[1] << endl;
				type1 = 0;  // adjusted for DISTANCE functions	
			}
			else if (type1 == 2) {  // query is RECTANGLE  
   			size1 = 4;			
				points1 = new double[size1];
				// Take Normalized Coords	
				points1[0] = ((Geo::Rectangle*) geo)->getP1()->getX(); 
   			points1[1] = ((Geo::Rectangle*) geo)->getP1()->getY(); 
   			points1[2] = ((Geo::Rectangle*) geo)->getP2()->getX(); 
   			points1[3] = ((Geo::Rectangle*) geo)->getP2()->getY();
				//cout << "nX_1: " << points1[0] << endl;  cout << "nY_1: " << points1[1] << endl;
				//cout << "nX_2: " << points1[2] << endl;  cout << "nY_2: " << points1[3] << endl;
			}
			// DEFINE 'QUERY' RELATED DATA < END >


			kNN_buffer.push_back(kNN_res());  // allocate a new entry to kNN buffer

			Result geoRes;
			geoRes.flags = Result::idAvailable;
			for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
				if ( (*iter)->subVar_exist ) 
					kNN_buffer[cnt].subject = (*iter)->value;  // get and store the subject id 		
				else if ( (*iter)->geoVar_exist )
					geoRes.id = (*iter)->value;  // get the ?geo unsigned id
				else 
					kNN_buffer[cnt].restArgs.push_back((*iter)->value);
			}

			kNN_buffer[cnt].geometry = geoRes.id;  // store the ?geo id 
			geoRes.ensureString(this);  // get the ?geo string value ('retrieve geometry')

			type2 = getLongLat(geoRes.value, points2, size2);  // set 'points2' & 'size2' values (NORMALIZED coords)

			// COMPUTE REAL DISTANCE < START >
			bool f = false;

			if (!type1 && !type2) {  // both are points
				//cout << "Point vs Point" << endl;
				kNN_buffer[cnt].distance = distance2(points1, points2);
			}
			else if ((f=(!type1 && type2==1)) || (!type2 && type1==1)) {  // point line
				if (f) {  // point line
					//cout << "Point vs Line" << endl;
					kNN_buffer[cnt].distance = pointLineStringDist2(points1, points2, size2);
				}
				else {  // line point
					//cout << "Line vs Point" << endl;
					kNN_buffer[cnt].distance = pointLineStringDist2(points2, points1, size1);
				}
			}
   		else if ((f=(!type1 && type2==2)) || (!type2 && type1==2)) {  // point polygon
				if (f) {  // point polygon
					//cout << "Point vs Polygon" << endl;
					kNN_buffer[cnt].distance = pointPolygonDist2(points1, points2, size2);
				}
				else {  // polygon point
					//cout << "Polygon vs Point" << endl;
					kNN_buffer[cnt].distance = pointPolygonDist2(points2, points1, size1);
				}	
   		}
			else if ((f=(type1==1 && type2==2)) || (type2==1 && type1==2)) {  // line polygon
				if (f) {  // line polygon
					//cout << "Line vs Polygon" << endl;
					kNN_buffer[cnt].distance = lineStringPolygonDist2(points1, size1, points2, size2);
				}
				else {  // polygon line
					//cout << "Polygon vs Line" << endl;
					kNN_buffer[cnt].distance = lineStringPolygonDist2(points2, size2, points1, size1);
				}
			}
			else if (type1==1 && type2==1) {  // line line
				//cout << "Line vs Line" << endl;	
				kNN_buffer[cnt].distance = lineStringLineStringDist2(points1, size1, points2, size2);
			}
			else if (type1==2 && type2==2) {  // polygon polygon
				//cout << "Polygon vs Polygon" << endl;
				kNN_buffer[cnt].distance = polygonPolygonDist2(points1, points2, size1, size2);
			}
			else if ((f=(!type1 && type2==3)) || (!type2 && type1==3)) {  // point multipoint
				if (f) {  // point multipoint	
					//cout << "Point vs Multipoint" << endl;
					kNN_buffer[cnt].distance = pointMultipointDist2(points1, points2, size2);
				}
				else {  // multipoint point
					//cout << "Multipoint vs Point" << endl;
					kNN_buffer[cnt].distance = pointMultipointDist2(points2, points1, size1);
				}
			}
			else if (type1==3 && type2==3) {  // multipoint multipoint
				//cout << "Multipoint vs Multipoint" << endl;
				kNN_buffer[cnt].distance = multipointMultipointDist2(points1, points2, size1, size2);	
			}
			else if ((f=(type1==1 && type2==3)) || (type2==1 && type1==3)) {  // line multipoint
				if (f) {  // line multipoint
					//cout << "Line vs Multipoint" << endl;
					kNN_buffer[cnt].distance = lineStringMultipointDist2(points1, size1, points2, size2);
				}
				else {  // multipoint line
					//cout << "Multipoint vs Line" << endl;
					kNN_buffer[cnt].distance = lineStringMultipointDist2(points2, size2, points1, size1);
				}
			}
			else if ((f=(type1==2 && type2==3)) || (type2==2 && type1==3)) {  // polygon multipoint
				if (f) {  // polygon multipoint
					//cout << "Polygon vs Multipoint" << endl;
					kNN_buffer[cnt].distance = polygonMultipointDist2(points1, points2, size1, size2);
				}
				else {  // multipoint polygon
					//cout << "Multipoint vs Polygon" << endl;
					kNN_buffer[cnt].distance = polygonMultipointDist2(points2, points1, size2, size1);	
				}
			}
			else {
				cout << "kNN Selection ERROR: Could not recognize pair of geometries: " << type1 << " " << type2 << endl;
				exit(-1);
			}
			// COMPUTE REAL DISTANCE < END >


			cnt++;

      	// Get the next tuples
			while (entry->next() != 0) {

				kNN_buffer.push_back(kNN_res());  // allocate a new entry to kNN buffer

				geoRes.flags = 0;  // reset the flags
				geoRes.flags = Result::idAvailable;
				for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
					if ( (*iter)->subVar_exist ) 
						kNN_buffer[cnt].subject = (*iter)->value;  // get and store the subject id 		
					else if ( (*iter)->geoVar_exist )
						geoRes.id = (*iter)->value;  // get the ?geo unsigned id	
					else 
						kNN_buffer[cnt].restArgs.push_back((*iter)->value);	
				}

				kNN_buffer[cnt].geometry = geoRes.id;  // store the ?geo id 
				geoRes.ensureString(this);  // get the ?geo string value ('retrieve geometry')		  
		
				type2 = getLongLat(geoRes.value, points2, size2);  // set 'points2' & 'size2' values (NORMALIZED coords)

				// COMPUTE REAL DISTANCE < START >
				f = false;

				if (!type1 && !type2) {  // both are points
					//cout << "Point vs Point" << endl;
					kNN_buffer[cnt].distance = distance2(points1, points2);
				}
				else if ((f=(!type1 && type2==1)) || (!type2 && type1==1)) {  // point line
					if (f) {  // point line
						//cout << "Point vs Line" << endl;
						kNN_buffer[cnt].distance = pointLineStringDist2(points1, points2, size2);
					}
					else {  // line point
						//cout << "Line vs Point" << endl;
						kNN_buffer[cnt].distance = pointLineStringDist2(points2, points1, size1);
					}
				}
   			else if ((f=(!type1 && type2==2)) || (!type2 && type1==2)) {  // point polygon
					if (f) {  // point polygon
						//cout << "Point vs Polygon" << endl;
						kNN_buffer[cnt].distance = pointPolygonDist2(points1, points2, size2);
					}
					else {  // polygon point
						//cout << "Polygon vs Point" << endl;
						kNN_buffer[cnt].distance = pointPolygonDist2(points2, points1, size1);
					}		
   			}
				else if ((f=(type1==1 && type2==2)) || (type2==1 && type1==2)) {  // line polygon
					if (f) {  // line polygon
						//cout << "Line vs Polygon" << endl;
						kNN_buffer[cnt].distance = lineStringPolygonDist2(points1, size1, points2, size2);
					}
					else {  // polygon line
						//cout << "Polygon vs Line" << endl;
						kNN_buffer[cnt].distance = lineStringPolygonDist2(points2, size2, points1, size1);
					}
				}
				else if (type1==1 && type2==1) {  // line line
					//cout << "Line vs Line" << endl;	
					kNN_buffer[cnt].distance = lineStringLineStringDist2(points1, size1, points2, size2);
				}
				else if (type1==2 && type2==2) {  // polygon polygon
					//cout << "Polygon vs Polygon" << endl;
					kNN_buffer[cnt].distance = polygonPolygonDist2(points1, points2, size1, size2);
				}
				else if ((f=(!type1 && type2==3)) || (!type2 && type1==3)) {  // point multipoint
					if (f) {  // point multipoint	
						//cout << "Point vs Multipoint" << endl;
						kNN_buffer[cnt].distance = pointMultipointDist2(points1, points2, size2);
					}
					else {  // multipoint point
						//cout << "Multipoint vs Point" << endl;
						kNN_buffer[cnt].distance = pointMultipointDist2(points2, points1, size1);
					}
				}
				else if (type1==3 && type2==3) {  // multipoint multipoint
					//cout << "Multipoint vs Multipoint" << endl;
					kNN_buffer[cnt].distance = multipointMultipointDist2(points1, points2, size1, size2);	
				}
				else if ((f=(type1==1 && type2==3)) || (type2==1 && type1==3)) {  // line multipoint
					if (f) {  // line multipoint
						//cout << "Line vs Multipoint" << endl;
						kNN_buffer[cnt].distance = lineStringMultipointDist2(points1, size1, points2, size2);
					}
					else {  // multipoint line
						//cout << "Multipoint vs Line" << endl;
						kNN_buffer[cnt].distance = lineStringMultipointDist2(points2, size2, points1, size1);
					}
				}
				else if ((f=(type1==2 && type2==3)) || (type2==2 && type1==3)) {  // polygon multipoint
					if (f) {  // polygon multipoint
						//cout << "Polygon vs Multipoint" << endl;
						kNN_buffer[cnt].distance = polygonMultipointDist2(points1, points2, size1, size2);
					}
					else {  // multipoint polygon
						//cout << "Multipoint vs Polygon" << endl;
						kNN_buffer[cnt].distance = polygonMultipointDist2(points2, points1, size2, size1);	
					}
				}	
				else {
					cout << "kNN Selection ERROR: Could not recognize pair of geometries: " << type1 << " " << type2 << endl;
					exit(-1);
				}
				// COMPUTE REAL DISTANCE < END >

	  			cnt++;	
			} 


			// Sort the buffer in ascending distance order
   		sort(kNN_buffer.begin(), kNN_buffer.end(), sortByDistance());
  

			// Testing Sort
			#if NEED_FOR_CODE
				cout << "Dist[0]: "  << fixed << kNN_buffer[0].distance << endl;
				cout << "Dist[1]: "  << fixed << kNN_buffer[1].distance << endl;
				cout << "Dist[2]: "  << fixed << kNN_buffer[2].distance << endl;
				cout << "Dist[3]: "  << fixed << kNN_buffer[3].distance << endl;
				cout << "Dist[4]: "  << fixed << kNN_buffer[4].distance << endl;	
				cout << "Dist[5]: "  << fixed << kNN_buffer[5].distance << endl;
				cout << "Dist[6]: "  << fixed << kNN_buffer[6].distance << endl;
				cout << "Dist[7]: "  << fixed << kNN_buffer[7].distance << endl;
				cout << "Dist[8]: "  << fixed << kNN_buffer[8].distance << endl;
				cout << "Dist[9]: "  << fixed << kNN_buffer[9].distance << endl;	
				cout << "Dist[10]: " << fixed << kNN_buffer[10].distance << endl;
				cout << "Dist[11]: " << fixed << kNN_buffer[11].distance << endl;
				cout << "Dist[12]: " << fixed << kNN_buffer[12].distance << endl;
				cout << "Dist[13]: " << fixed << kNN_buffer[13].distance << endl;
				cout << "Dist[14]: " << fixed << kNN_buffer[14].distance << endl;	
				cout << "Dist[15]: " << fixed << kNN_buffer[15].distance << endl;
				cout << "Dist[16]: " << fixed << kNN_buffer[16].distance << endl;
				cout << "Dist[17]: " << fixed << kNN_buffer[17].distance << endl;
				cout << "Dist[18]: " << fixed << kNN_buffer[18].distance << endl;
				cout << "Dist[19]: " << fixed << kNN_buffer[19].distance << endl;	
				//cout << "Taken Entities: " << cnt << endl;
			#endif

			// PASS REGISTER VALUES TO NEXT OPERATORS < START >
			cnt = 0;  pos = 0;
			for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
				if ((*iter)->select_exist) {
					if ( (*iter)->subVar_exist ) 
						(*iter)->value = kNN_buffer[cnt].subject;
					else if ( (*iter)->geoVar_exist )
						(*iter)->value = kNN_buffer[cnt].geometry;
					else  
						(*iter)->value = kNN_buffer[cnt].restArgs[pos++]; 
				}
			}
			// PASS REGISTER VALUES TO NEXT OPERATORS < END >

  			observedOutputCardinality++;
  			cnt++;

  			return 1;

  		#endif		
	}

#else

	#if RTREE_USED_IN_PLAN

		#if ONLY_POINT_GEOMS	

			//cout << "first...WITHIN RTREE Selection Evaluation (POINT_GEOMS)" << endl;

			observedOutputCardinality = 0;
			verif = filtered = 0;

			// Get the first tuple
			unsigned count = entry_RW->first();					
			if (!count) return 0;	

			observedOutputCardinality++;
			return 1;

		#else

			//cout << "first...WITHIN RTREE Selection Evaluation (ALL_GEOMS)" << endl;

			observedOutputCardinality = 0;
			verif = filtered = 0;

			// Get the first tuple
			unsigned count = entry_RW->first();
			if (!count) return 0;

			// DEFINE 'QUERY' RELATED DATA < START >
			type1_RW = geo_RW->type;
			if (type1_RW == 1) {  // query is POINT
				size1_RW = 2;
				points1_RW = new double[size1_RW];
				// Take Normalized Coords
				points1_RW[0] = ((Geo::Point*) geo_RW)->getX();
				points1_RW[1] = ((Geo::Point*) geo_RW)->getY();  	
			}
			else if (type1_RW == 2) {  // query is RECTANGLE  
   			size1_RW = 4;			
				points1_RW = new double[size1_RW];
				// Take Normalized Coords	
				points1_RW[0] = ((Geo::Rectangle*) geo_RW)->getP1()->getX(); 
   			points1_RW[1] = ((Geo::Rectangle*) geo_RW)->getP1()->getY(); 
   			points1_RW[2] = ((Geo::Rectangle*) geo_RW)->getP2()->getX(); 
   			points1_RW[3] = ((Geo::Rectangle*) geo_RW)->getP2()->getY();
			}
			// DEFINE 'QUERY' RELATED DATA < END >

			// Store the first tuple in priority queue 'SED'
			for (vector<Register*>::const_iterator iter = entryTail_RW.begin(), limit = entryTail_RW.end(); iter != limit; ++iter) {
				if ( (*iter)->subVar_exist ) 
					prRtreeEntry_RW.subject = (*iter)->value;  // get the subject id
				else if ( (*iter)->geoVar_exist )
					continue;
				else 
					prRtreeEntry_RW.restArgs.push_back((*iter)->value);
			}
			SED.push(prRtreeEntry_RW);

			// Get the ID of <hasGeometry> predicate
			#if LGD
				string hasGeom_string = "hasGeometry";
			#elif YAGO
				string hasGeom_string = "http://yago-knowledge.org/resource/hasGeometry";
			#endif
			Type::ID hasGeom_type;
			unsigned hasGeom_subType;
			bool exist;

			exist = runtime.getDatabase().getDictionary().lookup(hasGeom_string, hasGeom_type, hasGeom_subType, hasGeom_id_RW);
			if (exist) {
				//cout << "HASGEO_ID: " << hasGeom_id_RW << endl;
			}
			else { 
				cout << "HASGEO_ID NOT FOUND ERROR!" << endl;
				exit(-1);
			}

			// Get all the next tuples and store them in priority queue 'SED'
			while (entry_RW->next() != 0) {
				// Clear previous 'argument values'
				prRtreeEntry_RW.restArgs.clear();
				// Get the next tuples
				for (vector<Register*>::const_iterator iter = entryTail_RW.begin(), limit = entryTail_RW.end(); iter != limit; ++iter) {
					if ( (*iter)->subVar_exist ) 
						prRtreeEntry_RW.subject = (*iter)->value;  // get the subject id
					else if ( (*iter)->geoVar_exist )
						continue;
					else 
						prRtreeEntry_RW.restArgs.push_back((*iter)->value);
				}
				SED.push(prRtreeEntry_RW);	 
			}

			return next();

		#endif

	#else 

   	//cout << "first...Selection Evaluation" << endl;

   	observedOutputCardinality = 0;
   	predicate->setSelection(this);
   
   	sum = 0;

   	#if SEMIJOIN
		distThreshold = ((Selection::Distance*)predicate)->distThreshold;
   	#endif

   	verif = filtered = 0;

   	// Get the first tuple
   	unsigned count = input->first();
		if (!count) return 0;

   	// Match?
   	if (predicate->check()) {
    		observedOutputCardinality += count;
	  	 	sum += (count-1);
     	 	return count;
   	}

   	// Get the next one
   	return next();

   #endif	

#endif
}
//---------------------------------------------------------------------------
void Selection::getRegion(double *points, unsigned size, double &minX, double &maxX, double &minY, double &maxY)
   // Get the coordinates of the MBR for a given geometry
{
	minX = minY = 1000000;	//just a big value
	maxX = maxY = 0;	//the smallest possible, given that the coordinates are normalized in [0,360]

	//cout << "Size: " << size << endl;

	for(unsigned i=0;i<size;i++){

		if(i%2){//latitude

			if(minY>points[i])
				minY = points[i];

			if(maxY<points[i])
				maxY = points[i];
		}
		else{//longitude

			if(minX>points[i])
				minX = points[i];

			if(maxX<points[i])
				maxX = points[i];	
		}
	}

	//cout << "minX: " << minX << " maxX: " << maxX << " minY: " << minY << " maxY: " << maxY << endl;
	//getchar();
}
//---------------------------------------------------------------------------
unsigned Selection::getLongLat(const string &geometry, double* &points, unsigned &size) 
	// Return the type of geometry (0:point, 1:line, 2:polygon, 3:multipoint)
{
	size = (geometry.length() - sizeof(unsigned))/sizeof(double);
	if (size % 2) {
		cout << "Selection: Sth is wrong with the size of the geometry. Size is " << size << endl;
		exit(-1);
	}
	const char *data = geometry.data();
	points = (double*) &data[sizeof(unsigned)];
	return *(unsigned*) data;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#if K_NEAREST_NEIGHBOR

unsigned long long Selection::getHilbertId(double normX1, double normX2, double normY1, double normY2, double normCellSide, 
														 unsigned **c2i, unsigned cellsPerDim, unsigned b1, unsigned b2, unsigned &level)
	// Return the cell id
	// If the entity can be approximated only by the whole grid, then it returns the MAX_CELL_ID
	// The level of the cell is returned through the "level" argument
{
	unsigned prefX1, prefX2;  // minX, maxX
	unsigned prefY1, prefY2;  // minY, maxY

	if (normX1 == 1) prefX1 = cellsPerDim - 1;
	else prefX1 = normX1 / normCellSide;

	if (normX2 == 1) prefX2 = cellsPerDim - 1;
	else prefX2 = normX2 / normCellSide;

	if (normY1 == 1) prefY1 = cellsPerDim - 1;
	else prefY1 = normY1 / normCellSide;

	if (normY2 == 1) prefY2 = cellsPerDim - 1;
	else prefY2 = normY2 / normCellSide;

	// bottom-level cell
	if (prefX1 == prefX2 && prefY1 == prefY2) {

		level = b1 / 2;  // typically 13 for 32-bit integers

		// return the Hilbert-based cell id
		return c2i[prefX1][prefY1];		
	}
	else {  // upper-level cell

		//cout << "Computing the common father..." << endl;

		// find the common father of the lower left and the upper right cells
	
		// the cell ids
		unsigned ll = c2i[prefX1][prefY1],  // low-left bottom-level cell
			 ur = c2i[prefX2][prefY2],	      // upper-right bottom-level cell
			 mask = pow(2, b1-1), 		      // set left-most bit (cell IDs can have at most b1 bits)
			 bit1 = mask & ll, 		         // get the first bit of the low-left cell
			 bit2 = mask & ur,		         // get the first bit of the upper right cell
			 cnt = 0;

		//cout << "ll: " << ll << " ur: " << ur << endl;

		while (bit1 == bit2) {  // bits are the same

			// shift bit to the right
			mask >>= 1;

			// get the next bits
			bit1 = mask & ll;
			bit2 = mask & ur;

			// count common bits
			cnt++;
		}

		// handle odd numbers implicitly (the ID of the least common ancestor-cell is always defined by even number of bits)
		level = cnt / 2;

		if (!level) {  // the least common ancestor cell is the whole grid

			unsigned base = pow(2, b1), MAX_CELL_ID = base;
			
			while (base >= 4) {
				base >>= 2;
				MAX_CELL_ID += base;
			}
			MAX_CELL_ID--;

			// return the id corresponding to the whole grid
			return MAX_CELL_ID;
		}

		// return the upper-level cell id
		return (ll >> (b1 - 2*level));
	}
}
//---------------------------------------------------------------------------
unsigned Selection::getGridPos(unsigned cId, unsigned level, unsigned size, unsigned MAX_LEVEL, unsigned b1)
	// Return the "g_pos" of the cell with prefix "cId" at the given "level"
{
	//cout << "Getting grid pos." << endl;

	unsigned g_pos;

	if (!level) 
		g_pos = size - 1;  // max possible cell id (not to be confused with the max cell id at MAX_LEVEL)
	else {
		g_pos = cId;
		unsigned step = MAX_LEVEL;
		while (step > level) {
			g_pos += ((unsigned) 1 << (b1 - 2*(MAX_LEVEL - step)));  // pow(2, b1 - 2*(MAX_LEVEL - step));
			step--;
		} 
	}

	//cout << "End." << endl;

	return g_pos;
}
//---------------------------------------------------------------------------
bool Selection::getCellId(unsigned value, unsigned &c_id, unsigned &g_pos, unsigned &level)
	// Extract the cell id "c_id", the grid pos "g_pos" and the level "level" from the ID of an entity 
	// Return "true" if the input entity is associated with spatial info, "false" otherwise
{
	level = 0;
	c_id = value;

	if (c_id & 1) 
	   level = ((c_id & LVMASK) >> 1);
	else 
      return false;
	
	c_id >>= (grid->b - 2*level);

	if (!level) 
		c_id = g_pos = (grid->size - 1);
	else {
		g_pos = c_id;
		unsigned step = grid->MAX_LEVEL;
		while (step > level) {
			g_pos += ((unsigned) 1 << (grid->b1 - 2*(grid->MAX_LEVEL - step)));  // pow(2, b1 - 2*(MAX_LEVEL - step));
			step--;
		} 
	}

	return true;
}
//---------------------------------------------------------------------------
void Selection::parseGeometry(Geo::Geometry *g, double* &points, unsigned &size)
	// Get the geometry and assign points and size values 
{
	if (g->type == 1) {  // Geometry: POINT
		size = 2;
		points = new double[size];

		stringstream ss1(stringstream::in | stringstream::out), ss2(stringstream::in | stringstream::out);

		ss1 << dynamic_cast<Geo::Point*>(g)->getX();
		ss1 >> points[0];
		cout << "POINT, getX() is " << points[0] << endl;  
		ss2 << dynamic_cast<Geo::Point*>(g)->getY();
		ss2 >> points[1];
		cout << "POINT, getY() is " << points[1] << endl;
	}
	else if (g->type == 2) {  // Geometry: RECTANGLE
		size = 4;	
		points = new double[size];

		stringstream ss1(stringstream::in | stringstream::out), ss2(stringstream::in | stringstream::out),
						 ss3(stringstream::in | stringstream::out), ss4(stringstream::in | stringstream::out);
 	
		ss1 << dynamic_cast<Geo::Rectangle*>(g)->getP1()->getX();
		ss1 >> points[0];
		cout << "RECTANGLE, getP1()->getX() is " << points[0] << endl;
		ss2 << dynamic_cast<Geo::Rectangle*>(g)->getP1()->getY();
		ss2 >> points[1];
		cout << "RECTANGLE, getP1()->getY() is " << points[1] << endl;

		ss3 << dynamic_cast<Geo::Rectangle*>(g)->getP2()->getX();
		ss3 >> points[2];
		cout << "RECTANGLE, getP2()->getX() is " << points[2] << endl;
		ss4 << dynamic_cast<Geo::Rectangle*>(g)->getP2()->getY();
		ss4 >> points[3];
		cout << "RECTANGLE, getP2()->getY() is " << points[3] << endl;
	}
}
//---------------------------------------------------------------------------
double Selection::RtreeDistPointRegion(double *point, double *region, unsigned size)
	// Compute the square of the minimum distance between a point and a region,
	// based on the Rtree Point-Region minimum distance libspatialindex library
{
	double ret = 0.0;

	for (unsigned i = 0; i < size; i++)
	{
		if (point[i] < region[i])
		{
			ret += pow(region[i] - point[i], 2.0);
		}
		else if (point[i] > region[i+2])
		{
			ret += pow(point[i] - region[i+2], 2.0);
		}
	}

	return ret;
}
//---------------------------------------------------------------------------	
void Selection::findEntityLimits(unsigned *limits, unsigned* &grid, unsigned maxCell, unsigned b2)
	// Find the 'limit entity conditions' for each level based on 'maxCell'
{
	unsigned prevCells;

	// L_13 limit
	prevCells = 0;
	limits[0] = 0 << 31;   
	if (grid[maxCell] != 0) {
		limits[0] |= 1;
		limits[0] |= (13 << 1); 	
		limits[0] |= ((grid[maxCell] - 1) << (b2 - 1));
		limits[0] |= (maxCell << b2);
	}
	else {  // create an artificial spatial entity
		limits[0] |= 1;
		limits[0] |= (13 << 1);	
		limits[0] |= (0 << (b2 - 1));		
		limits[0] |= (maxCell << b2);		
	}

	// L_12 limit
	prevCells += 8192*8192;
	limits[1] = 0 << 31;	 	
	if (grid[prevCells + (maxCell >> 2)] != 0) {	
		limits[1] |= 1;
		limits[1] |= (12 << 1);
		limits[1] |= ((grid[prevCells + (maxCell >> 2)] - 1) << (b2 - 1)); 
		limits[1] |= ((maxCell >> 2) << (b2 + 2)); 
	}
	else {  // create an artificial spatial entity
		limits[1] |= 1;
		limits[1] |= (12 << 1);	
		limits[1] |= (0 << (b2 - 1));
		limits[1] |= ((maxCell >> 2) << (b2 + 2));
	}

	// L_11 limit	
	prevCells += 4096*4096;
	limits[2] = 0 << 31;	 
	if (grid[prevCells + (maxCell >> 4)] != 0) {	
		limits[2] |= 1;
		limits[2] |= (11 << 1);
		limits[2] |= ((grid[prevCells + (maxCell >> 4)] - 1) << (b2 - 1));
		limits[2] |= ((maxCell >> 4) << (b2 + 4));
	}
	else {  // create an artificial spatial entity
		limits[2] |= 1;
		limits[2] |= (11 << 1);
		limits[2] |= (0 << (b2 - 1));
		limits[2] |= ((maxCell >> 4) << (b2 + 4));
	}

	// L_10 limit 	
	prevCells += 2048*2048;
	limits[3] = 0 << 31;	 
	if (grid[prevCells + (maxCell >> 6)] != 0) {	
		limits[3] |= 1;
		limits[3] |= (10 << 1);		
		limits[3] |= ((grid[prevCells + (maxCell >> 6)] - 1) << (b2 - 1));
		limits[3] |= ((maxCell >> 6) << (b2 + 6));
	}
	else {  // create an artificial spatial entity
		limits[3] |= 1;
		limits[3] |= (10 << 1);	
		limits[3] |= (0 << (b2 - 1));
		limits[3] |= ((maxCell >> 6) << (b2 + 6));
	}

	// L_9 limit
	prevCells += 1024*1024;
	limits[4] = 0 << 31;
	if (grid[prevCells + (maxCell >> 8)] != 0) {
		limits[4] |= 1;
		limits[4] |= (9 << 1);
		limits[4] |= ((grid[prevCells + (maxCell >> 8)] - 1) << (b2 - 1));
		limits[4] |= ((maxCell >> 8) << (b2 + 8));
	}
	else {  // create an artificial spatial entity
		limits[4] |= 1;
		limits[4] |= (9 << 1);	
		limits[4] |= (0 << (b2 - 1));
		limits[4] |= ((maxCell >> 8) << (b2 + 8));
	}

	// L_8 limit
	prevCells += 512*512;
	limits[5] = 0 << 31;
	if (grid[prevCells + (maxCell >> 10)] != 0) {
		limits[5] |= 1;
		limits[5] |= (8 << 1);	
		limits[5] |= ((grid[prevCells + (maxCell >> 10)] - 1) << (b2 - 1));	
		limits[5] |= ((maxCell >> 10) << (b2 + 10));
	}
	else {  // create an artificial spatial entity
		limits[5] |= 1;
		limits[5] |= (8 << 1);
		limits[5] |= (0 << (b2 - 1));
		limits[5] |= ((maxCell >> 10) << (b2 + 10));
	}

	// L_7 limit
	prevCells += 256*256;
	limits[6] = 0 << 31;	
	if (grid[prevCells + (maxCell >> 12)] != 0) {
		limits[6] |= 1;
		limits[6] |= (7 << 1);
		limits[6] |= ((grid[prevCells + (maxCell >> 12)] - 1) << (b2 - 1));	
		limits[6] |= ((maxCell >> 12) << (b2 + 12));	
	}
	else {  // create an artificial spatial entity
		limits[6] |= 1;
		limits[6] |= (7 << 1);
		limits[6] |= (0 << (b2 - 1));
		limits[6] |= ((maxCell >> 12) << (b2 + 12));
	}	
	
	// L_6 limit
	prevCells += 128*128;
	limits[7] = 0 << 31;
	if (grid[prevCells + (maxCell >> 14)] != 0) {
		limits[7] |= 1;
		limits[7] |= (6 << 1);	
		limits[7] |= ((grid[prevCells + (maxCell >> 14)] - 1) << (b2 - 1));	
		limits[7] |= ((maxCell >> 14) << (b2 + 14));
	}
	else {  // create an artificial spatial entity
		limits[7] |= 1;
		limits[7] |= (6 << 1);
		limits[7] |= (0 << (b2 - 1));
		limits[7] |= ((maxCell >> 14) << (b2 + 14));
	}	

	// L_5 limit	
	prevCells += 64*64;
	limits[8] = 0 << 31;	
	if (grid[prevCells + (maxCell >> 16)] != 0) {
		limits[8] |= 1;
		limits[8] |= (5 << 1);
		limits[8] |= ((grid[prevCells + (maxCell >> 16)] - 1) << (b2 - 1));	
		limits[8] |= ((maxCell >> 16) << (b2 + 16));	
	}	
	else {  // create an artificial spatial entity
		limits[8] |= 1;
		limits[8] |= (5 << 1);
		limits[8] |= (0 << (b2 - 1));
		limits[8] |= ((maxCell >> 16) << (b2 + 16));	
	}

	// L_4 limit
	prevCells += 32*32;
	limits[9] = 0 << 31;	
	if (grid[prevCells + (maxCell >> 18)] != 0) {
		limits[9] |= 1;
		limits[9] |= (4 << 1);
		limits[9] |= ((grid[prevCells + (maxCell >> 18)] - 1) << (b2 - 1));	
		limits[9] |= ((maxCell >> 18) << (b2 + 18));
	}	
	else {  // create an artificial spatial entity
		limits[9] |= 1;
		limits[9] |= (4 << 1);
		limits[9] |= (0 << (b2 - 1));
		limits[9] |= ((maxCell >> 18) << (b2 + 18));
	}

	// L_3 limit
	prevCells += 16*16;
	limits[10] = 0 << 31;	
	if (grid[prevCells + (maxCell >> 20)] != 0) {
		limits[10] |= 1;
		limits[10] |= (3 << 1);	
		limits[10] |= ((grid[prevCells + (maxCell >> 20)] - 1) << (b2 - 1));
		limits[10] |= ((maxCell >> 20) << (b2 + 20));	
	}
	else {  // create an artificial spatial entity
		limits[10] |= 1;
		limits[10] |= (3 << 1);
		limits[10] |= (0 << (b2 - 1));
		limits[10] |= ((maxCell >> 20) << (b2 + 20));	
	}

	// L_2 limit
	prevCells += 8*8;
	limits[11] = 0 << 31;
	if (grid[prevCells + (maxCell >> 22)] != 0) {
		limits[11] |= 1;
		limits[11] |= (2 << 1);	
		limits[11] |= ((grid[prevCells + (maxCell >> 22)] - 1) << (b2 - 1));	
		limits[11] |= ((maxCell >> 22) << (b2 + 22));	
	}
	else {  // create an artificial spatial entity
		limits[11] |= 1;
		limits[11] |= (2 << 1);	
		limits[11] |= (0 << (b2 - 1));
		limits[11] |= ((maxCell >> 22) << (b2 + 22));
	}

	// L_1 limit
	prevCells += 4*4;	
	limits[12] = 0 << 31;
	if (grid[prevCells + (maxCell >> 24)] != 0) {	
		limits[12] |= 1;		
		limits[12] |= (1 << 1);
		limits[12] |= ((grid[prevCells + (maxCell >> 24)] - 1) << (b2 - 1));	
		limits[12] |= ((maxCell >> 24) << (b2 + 24));
	}
	else {  // create an artificial spatial entity
		limits[12] |= 1;
		limits[12] |= (1 << 1);
		limits[12] |= (0 << (b2 - 1));
		limits[12] |= ((maxCell >> 24) << (b2 + 24));
	}
}
//---------------------------------------------------------------------------
unsigned Selection::findMaxCellInRecZone(unsigned zone, unsigned dir, unsigned** &c2i, double *q_ll, double cSide) 
	// Find the maximum bottom cell id of '(zone, dir)' rectangle zone
{
	double DIRECT_ll[2];
	unsigned DIRECT_cell;
	unsigned max_cell, next_cell, prev_cell;
	double TEMP_ll[2];

	// Find the 'direct' direction 'low-left' values & 'cells' (based on 'zone')
	if (dir == 0) {  // U  
		DIRECT_ll[0] = q_ll[0];
		DIRECT_ll[1] = q_ll[1] + (zone * cSide);	

		if ((DIRECT_ll[1] + cSide) > 1.0)	
			DIRECT_cell = gridCells;  // out-of-grid value
		else 
			DIRECT_cell = c2i[(unsigned)(DIRECT_ll[0] / cSide)][(unsigned)(DIRECT_ll[1] / cSide)];
		//cout << "U_cell: " << DIRECT_cell << endl;	
	}

	else if (dir == 1) {  // R 
		DIRECT_ll[0] = q_ll[0] + (zone * cSide);
		DIRECT_ll[1] = q_ll[1]; 	

		if ((DIRECT_ll[0] + cSide) > 1.0)
			DIRECT_cell = gridCells;  // out-of-grid value
		else 	
			DIRECT_cell = c2i[(unsigned)(DIRECT_ll[0] / cSide)][(unsigned)(DIRECT_ll[1] / cSide)];
		//cout << "R_cell: " << DIRECT_cell << endl;
	}

	else if (dir == 2) {  // D
		DIRECT_ll[0] = q_ll[0];
		DIRECT_ll[1] = q_ll[1] - (zone * cSide);		

		if (DIRECT_ll[1] < 0.0)
			DIRECT_cell = gridCells;  // out-of-grid value
		else
			DIRECT_cell = c2i[(unsigned)(DIRECT_ll[0] / cSide)][(unsigned)(DIRECT_ll[1] / cSide)];
		//cout << "D_cell: " << DIRECT_cell << endl;
	}

	else if (dir == 3) {  // L
		DIRECT_ll[0] = q_ll[0] - (zone * cSide);
		DIRECT_ll[1] = q_ll[1];

		if (DIRECT_ll[0] < 0.0)	
			DIRECT_cell = gridCells;  // out-of-grid value
		else 
			DIRECT_cell = c2i[(unsigned)(DIRECT_ll[0] / cSide)][(unsigned)(DIRECT_ll[1] / cSide)];
		//cout << "L_cell: " << DIRECT_cell << endl;
	} 


	if (DIRECT_cell == gridCells) {  // Rectangle Zone is Out-Of-Grid 
		max_cell = gridCells; 
		return max_cell;
	}
	else {  // Rectangle Zone is Within-Grid
		max_cell = DIRECT_cell;
		// the process continues below... 
	}
	

	// FIND THE 'max_cell' OF RESPECTIVE RECTANGLE ZONE < START > 
	if (dir == 0) {  // U
		// RIGHT direction
		next_cell = DIRECT_cell;
		TEMP_ll[0] = DIRECT_ll[0] + cSide;
		TEMP_ll[1] = DIRECT_ll[1];
		while (1) {
			// Get the 'next RIGHT cell'
			if ( (TEMP_ll[0] + cSide) > 1.0 || (TEMP_ll[0] - DIRECT_ll[0]) > (zone * cSide) ) 
				break;  // 'out-of-grid' OR 'out-of-zone' cell	
			else
				next_cell = c2i[(unsigned)(TEMP_ll[0] / cSide)][(unsigned)(TEMP_ll[1] / cSide)];
			//cout << "U_next_cell: " << next_cell << endl;
			// Update for next iteration
			TEMP_ll[0] = TEMP_ll[0] + cSide;
			// Update max value
			if (next_cell > max_cell) 
				max_cell = next_cell;
		}
		// LEFT direction
		prev_cell = DIRECT_cell;
		TEMP_ll[0] = DIRECT_ll[0] - cSide;
		TEMP_ll[1] = DIRECT_ll[1];	
		while (1) {
			// Get the 'next LEFT cell'	
			if ( TEMP_ll[0] < 0.0 || (DIRECT_ll[0] - TEMP_ll[0]) > ((zone - 1) * cSide) ) 
				break;  // 'out-of-grid' OR 'out-of-zone' cell
			else
				prev_cell = c2i[(unsigned)(TEMP_ll[0] / cSide)][(unsigned)(TEMP_ll[1] / cSide)];
			//cout << "U_prev_cell: " << prev_cell << endl;
			// Update for next iteration
			TEMP_ll[0] = TEMP_ll[0] - cSide;
			// Update max value
			if (prev_cell > max_cell)	
				max_cell = prev_cell;
		}
		//cout << endl << "U_max_cell: " << max_cell << endl << endl;
	}

	else if (dir == 1) {  // R 	
		// UP direction
		next_cell = DIRECT_cell; 
		TEMP_ll[0] = DIRECT_ll[0];
		TEMP_ll[1] = DIRECT_ll[1] + cSide;
		while (1) {
			// Get the 'next UP cell' 		
			if ( (TEMP_ll[1] + cSide) > 1.0 || (TEMP_ll[1] - DIRECT_ll[1]) > ((zone - 1) * cSide) )			
				break;  // 'out-of-grid' OR 'out-of-zone' cell
			else		
				next_cell = c2i[(unsigned)(TEMP_ll[0] / cSide)][(unsigned)(TEMP_ll[1] / cSide)];
			//cout << "R_next_cell: " << next_cell << endl;
			// Update for next iteration
			TEMP_ll[1] = TEMP_ll[1] + cSide;	
			// Update max value	
			if (next_cell > max_cell)
				max_cell = next_cell;
		}
		// DOWN direction
		prev_cell = DIRECT_cell;
		TEMP_ll[0] = DIRECT_ll[0];
		TEMP_ll[1] = DIRECT_ll[1] - cSide;
		while (1) {
			// Get the 'next DOWN cell'	 		
			if ( TEMP_ll[1] < 0.0 || (DIRECT_ll[1] - TEMP_ll[1]) > (zone * cSide) )			
				break;  // 'out-of-grid' OR 'out-of-zone' cell
			else		
				prev_cell = c2i[(unsigned)(TEMP_ll[0] / cSide)][(unsigned)(TEMP_ll[1] / cSide)];
			//cout << "R_prev_cell: " << prev_cell << endl;
			// Update for next iteration
			TEMP_ll[1] = TEMP_ll[1] - cSide;
			// Update max value
			if (prev_cell > max_cell)
				max_cell = prev_cell;
		}	
		//cout << "R_max_cell: " << max_cell << endl << endl;
	}

	else if (dir == 2) {  // D
		// RIGHT direction
		next_cell = DIRECT_cell;	
		TEMP_ll[0] = DIRECT_ll[0] + cSide;
		TEMP_ll[1] = DIRECT_ll[1];
		while (1) {
			// Get the 'next RIGHT cell'
			if ( (TEMP_ll[0] + cSide) > 1.0 || (TEMP_ll[0] - DIRECT_ll[0]) > ((zone - 1) * cSide) )			
				break;  // 'out-of-grid' OR 'out-of-zone' cell
			else		
				next_cell = c2i[(unsigned)(TEMP_ll[0] / cSide)][(unsigned)(TEMP_ll[1] / cSide)];
			//cout << "D_next_cell: " << next_cell << endl;
			// Update for next iteration
			TEMP_ll[0] = TEMP_ll[0] + cSide;		
			// Update max value
			if (next_cell > max_cell)
				max_cell = next_cell;
		}
		// LEFT direction
		prev_cell = DIRECT_cell;
		TEMP_ll[0] = DIRECT_ll[0] - cSide;
		TEMP_ll[1] = DIRECT_ll[1];
		while (1) {
			// Get the 'next LEFT cell'	
			if ( TEMP_ll[0] < 0.0 || (DIRECT_ll[0] - TEMP_ll[0]) > (zone * cSide) )
				break;  // 'out-of-grid' OR 'out-of-zone' cell	
			else
				prev_cell = c2i[(unsigned)(TEMP_ll[0] / cSide)][(unsigned)(TEMP_ll[1] / cSide)];
			//cout << "D_prev_cell: " << prev_cell << endl;
			// Update for next iteration
			TEMP_ll[0] = TEMP_ll[0] - cSide;
			// Update max value	
			if (prev_cell > max_cell)	
				max_cell = prev_cell;
		}
		//cout << "D_max_cell: " << max_cell << endl << endl;
	}

	else if (dir == 3) {  // L
		// UP direction
		next_cell = DIRECT_cell; 
		TEMP_ll[0] = DIRECT_ll[0];
		TEMP_ll[1] = DIRECT_ll[1] + cSide;
		while (1) {
			// Get the 'next UP cell' 		
			if ( (TEMP_ll[1] + cSide) > 1.0 || (TEMP_ll[1] - DIRECT_ll[1]) > (zone * cSide) )			
				break;  // 'out-of-grid' OR 'out-of-zone' cell
			else		
				next_cell = c2i[(unsigned)(TEMP_ll[0] / cSide)][(unsigned)(TEMP_ll[1] / cSide)];
			//cout << "L_next_cell: " << next_cell << endl;
			// Update for next iteration
			TEMP_ll[1] = TEMP_ll[1] + cSide;
			// Update max value	
			if (next_cell > max_cell)
				max_cell = next_cell;
		}
		// DOWN direction
		prev_cell = DIRECT_cell;
		TEMP_ll[0] = DIRECT_ll[0];
		TEMP_ll[1] = DIRECT_ll[1] - cSide;	
		while (1) {
			// Get the 'next DOWN cell'	
			if ( TEMP_ll[1] < 0.0 || (DIRECT_ll[1] - TEMP_ll[1]) > ((zone - 1) * cSide) )			
				break;  // 'out-of-grid' OR 'out-of-zone' cell
			else		
				prev_cell = c2i[(unsigned)(TEMP_ll[0] / cSide)][(unsigned)(TEMP_ll[1] / cSide)];
			//cout << "L_prev_cell: " << prev_cell << endl;
			// Update for next iteration	
			TEMP_ll[1] = TEMP_ll[1] - cSide;
			// Update max value
			if (prev_cell > max_cell)
				max_cell = prev_cell;	
		}
		//cout << "L_max_cell: " << max_cell << endl << endl;
	}
	// FIND THE 'max_cell' OF RESPECTIVE RECTANGLE ZONE < END >

	return max_cell;	
}
//---------------------------------------------------------------------------	
double Selection::dot(double *A, double *B)
	// Compute the dot product of vectors A and B
{
        return A[0] * B[0] + A[1] * B[1];
}
//---------------------------------------------------------------------------
double Selection::distance2(double *A, double *B)
	// Compute the square of the distance from A to B
{
        double d1 = A[0] - B[0];
        double d2 = A[1] - B[1];
        return d1*d1+d2*d2;
}
//---------------------------------------------------------------------------
double Selection::pointLineDist2(double *point, double *segment)
	// Compute the square of the minimum distance between a point and a line
{       
	double v[2],w[2];

	//Vector v = S.P1 - S.P0;
	v[0] = segment[2] - segment[0];
	v[1] = segment[3] - segment[1];

	//Vector w = P - S.P0;
	w[0] = point[0] - segment[0];
	w[1] = point[1] - segment[1];

	double c1 = dot(w,v);
     	if ( c1 <= 0 )
          	return distance2(point,segment);	//distance2(P, S.P0);

     	double c2 = dot(v,v);
     	if ( c2 <= c1 )
          	return distance2(point,&segment[2]);	//distance2(P, S.P1);

     	double b = c1 / c2;
	double pb[2];

	//Point Pb = S.P0 + b * v;
	pb[0] = segment[0] + b*v[0];
	pb[1] = segment[1] + b*v[1];

    	return distance2(point, pb);	//return distance2(P, Pb);
}
//---------------------------------------------------------------------------
double Selection::pointLineStringDist2(double *point, double *line, unsigned size)
	// Compute the square of the minimum distance between a point and a line segment
{       
	double d, mD = 100000000;	//the minimum distance
	
	size -= 2;

	//for each segment of the given linestring
	for(unsigned i=0;i<size;i+=2){

		//distance between the line and the point
		d = pointLineDist2(point,&line[i]);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::lineLineDist2(double *seg1, double *seg2)
	// Compute the square of the minimum distance between a line and a line
{ 
	double u[2],v[2],w[2];

	//Vector   u = S1.P1 - S1.P0;
	u[0] = seg1[2] - seg1[0];
	u[1] = seg1[3] - seg1[1];

	//Vector   v = S2.P1 - S2.P0;
	v[0] = seg2[2] - seg2[0];
	v[1] = seg2[3] - seg2[1];
	
	//Vector   w = S1.P0 - S2.P0;
	w[0] = seg1[0] - seg2[0];
	w[1] = seg1[1] - seg2[1];

	double    a = dot(u,u);         // always >= 0
	double    b = dot(u,v);
	double    c = dot(v,v);         // always >= 0
	double    d = dot(u,w);
	double    e = dot(v,w);
	double    D = a*c - b*b;        // always >= 0
	double    sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
	double    tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

	// compute the line parameters of the two closest points
	if (D < 0.00000001) { 		// the lines are almost parallel
		sN = 0.0;         	// force using point P0 on segment S1
		sD = 1.0;         	// to prevent possible division by 0.0 later
		tN = e;
		tD = c;
	}
	else {                 		// get the closest points on the infinite lines
		sN = (b*e - c*d);
		tN = (a*e - b*d);
		
		if (sN < 0.0) {        	// sc < 0 => the s=0 edge is visible
	    		sN = 0.0;
	    		tN = e;
	    		tD = c;
		}
		else if (sN > sD) {	// sc > 1  => the s=1 edge is visible
	    		sN = sD;
	    		tN = e + b;
	    		tD = c;
		}
	}

	if (tN < 0.0) {           	// tc < 0 => the t=0 edge is visible

		tN = 0.0;

		// recompute sc for this edge
		if (-d < 0.0)
	    		sN = 0.0;
		else if (-d > a)
	    		sN = sD;
		else {
	    		sN = -d;
	    		sD = a;
		}
	}
	else if (tN > tD) {      	// tc > 1  => the t=1 edge is visible

		tN = tD;
	
		// recompute sc for this edge
		if ((-d + b) < 0.0)
	    		sN = 0;
		else if ((-d + b) > a)
		    sN = sD;
		else {
	    		sN = (-d +  b);
	    		sD = a;
		}
	}

	// finally do the division to get sc and tc
	sc = (abs(sN) < 0.00000001 ? 0.0 : sN / sD);
	tc = (abs(tN) < 0.00000001 ? 0.0 : tN / tD);

	// get the difference of the two closest points
	double dp[2];

	//Vector   dP = w + (sc * u) - (tc * v);  // =  S1(sc) - S2(tc)
	dp[0] = w[0] + sc *u[0] - tc*v[0];
	dp[1] = w[1] + sc *u[1] - tc*v[1];

	return dot(dp,dp);   // return the square of the minimum distance
}
//---------------------------------------------------------------------------
double Selection::lineStringLineStringDist2(double *line1, unsigned size1, double *line2, unsigned size2)
	// Compute the square of the minimum distance between a line segment and a line segment
{       
	double d, mD = 100000000;	//the minimum distance
	
	size1 -= 2;
	size2 -= 2;

	//for each segment of the first linestring
	for(unsigned i=0;i<size1;i+=2){

		//for each segment of the second linestring
		for(unsigned j=0;j<size2;j+=2){

			//distance between the lines
			d = lineLineDist2(&line1[i],&line2[j]);

			if(!d) return 0;

			if(d<mD) mD = d;
		}
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::linePolygonDist2(double *line, double *polygon, unsigned size)
	// Compute the square of the minimum distance between a line and a polygon
{ 
	double d, mD = 100000000;	//the minimum distance
	
	size -= 2;

	//for each line of the given polygon
	for(unsigned i=0;i<size;i+=2){

		//distance between the lines
		d = lineLineDist2(line,&polygon[i]);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::lineStringPolygonDist2(double *line, unsigned size1, double *polygon, unsigned size2)
	// Compute the square of the minimum distance between a line segment and a polygon
{       
	double d, mD = 100000000;	//the minimum distance
	
	size1 -= 2;

	//for each segment of the given linestring
	for(unsigned i=0;i<size1;i+=2){

		d = linePolygonDist2(&line[i],polygon,size2);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::pointPolygonDist2(double *point, double *polygon, unsigned size)
	// Compute the square of the minimum distance between a point and a polygon
{ 
	double d, mD = 100000000;	//the minimum distance
	
	size -= 2;

	//for each line of the given polygon
	for(unsigned i=0;i<size;i+=2){

		//distance between the point and the line
		d = pointLineDist2(point,&polygon[i]);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::polygonPolygonDist2(double *polygon1, double *polygon2, unsigned size1, unsigned size2)
	// Compute the square of the minimum distance between a polygon and a polygon
{ 
	double d, mD = 100000000;	//the minimum distance
	
	size1 -= 2;
	size2 -= 2;

	//for each line of the first polygon
	for(unsigned i=0;i<size1;i+=2){

		//for each line of the second polygon
		for(unsigned j=0;j<size2;j+=2){

			//distance between the lines
			d = lineLineDist2(&polygon1[i],&polygon2[j]);

			if(!d) return 0;

			if(d<mD) mD = d;
		}
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::pointMultipointDist2(double *point, double *multipoint, unsigned size)
	// Compute the square of the minimum distance between a point and a multipoint
{ 
	double d, mD = 100000000;	//the minimum distance
	
	//for each point of the given multipoint
	for(unsigned i=0;i<size;i+=2){

		//distance between the points
		d = distance2(point,&multipoint[i]);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::lineMultipointDist2(double *line, double *multipoint, unsigned size)
    	//Compute the square of the minimum distance between a line and a multipoint
{ 
	double d, mD = 100000000;	//the minimum distance
	
	//for each point of the given multipoint
	for(unsigned i=0;i<size;i+=2){

		//distance between the point and the the line
		d = pointLineDist2(&multipoint[i],line);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::lineStringMultipointDist2(double *line, unsigned size1, double *multipoint, unsigned size2)
	// Compute the square of the minimum distance between a line segment and a multipoint
{ 
	double d, mD = 100000000;	//the minimum distance
	
	size1 -= 2;

	//for each segment of the given linestring
	for(unsigned i=0;i<size1;i+=2){

		//distance between the point and the the line
		d = lineMultipointDist2(&line[i],multipoint,size2);

		if(!d) return 0;

		if(d<mD) mD = d;
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::polygonMultipointDist2(double *polygon, double *multipoint, unsigned size1, unsigned size2)
	// Compute the square of the minimum distance between a polygon and a multipoint
{ 
	double d, mD = 100000000;	//the minimum distance
	
	size1 -= 2;

	//for each line of the given polygon
	for(unsigned i=0;i<size1;i+=2){

		//for each point of the given multipoint
		for(unsigned j=0;j<size2;j+=2){

			//distance between the point and the the line
			d = pointLineDist2(&multipoint[j],&polygon[i]);

			if(!d) return 0;

			if(d<mD) mD = d;
		}
	}

	return mD;
}
//---------------------------------------------------------------------------
double Selection::multipointMultipointDist2(double *multipoint1, double *multipoint2, unsigned size1, unsigned size2)
	// Compute the square of the minimum distance between a multipoint and a multipoint
{ 
	double d, mD = 100000000;	//the minimum distance
	
	//for each point of the first multipoint
	for(unsigned i=0;i<size1;i+=2){

		//for each point of the second multipoint
		for(unsigned j=0;j<size2;j+=2){

			//distance between the points
			d = distance2(&multipoint1[i],&multipoint2[j]);

			if(!d) return 0;

			if(d<mD) mD = d;
		}
	}

	return mD;
}

#endif

//---------------------------------------------------------------------------
#if SEMIJOIN
void Selection::lookUpGeo(unsigned id, string& g)
{
	    const char* start,*stop;
		/// The type
		Type::ID type;
		/// The sub-type
		unsigned subType;

	    if((itter=cache.find(id))!=cache.end()){
			g = itter->second;
	    }
	    else if(runtime.getDatabase().getDictionary().lookupById(id,start,stop,type,subType)){
		    g=string(start,stop);
		    cache.insert(pair<unsigned,string> (id,g));
	    }
}
//---------------------------------------------------------------------------
void Selection::getLongLat(string& geometry, double& ln, double& lt)
{
	//cout << "Geometry: " << geometry << endl;

	string dlm = ",";
	std::size_t pos = geometry.find(dlm);

	if (pos==std::string::npos){

		cout << "Comma NOT found in geometry!" << endl;
		exit(-1);
	}
	
	string lon = geometry.substr (0,pos);
	string lat = geometry.substr (pos+1);
	//cout << "Long: " << lon << endl;
	//cout << "Lat: " << lat << endl;

	stringstream ss1(stringstream::in|stringstream::out), ss2(stringstream::in|stringstream::out);
	ss1 << lon;
	ss1 >> ln;
	ss2 << lat;
	ss2 >> lt;
}
#endif

//---------------------------------------------------------------------------
unsigned Selection::next()
   // Produce the next tuple
{
#if K_NEAREST_NEIGHBOR

	if (!kNN_encode) 
	{
		//cout << "next...kNN BASELINE or BASIC Selection Evaluation" << endl;

   	// Fetch all the 'argument values' that exist in 'SELECT' query clause
   	if ( (cnt == k_nearest_neighbor) || (cnt == kNN_buffer.size()) )
      	return 0;      

		pos = 0;
		for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
			if ((*iter)->select_exist) {
				if ( (*iter)->subVar_exist ) 
					(*iter)->value = kNN_buffer[cnt].subject;
				else if ( (*iter)->geoVar_exist )
					(*iter)->value = kNN_buffer[cnt].geometry;
				else  
					(*iter)->value = kNN_buffer[cnt].restArgs[pos++]; 
			}
		}

   	observedOutputCardinality++;
   	cnt++;

   	return 1;
	}

	else 
	{
		//cout << "next...kNN ENCODING Selection Evaluation" << endl;	
		
   	// Fetch all the 'argument values' that exist in 'SELECT' query clause
   	if (cnt == kNN_buffer.size())
      	return 0; 

		pos = 0;
		for (vector<Register*>::const_iterator iter = entryTail.begin(), limit = entryTail.end(); iter != limit; ++iter) {
			if ((*iter)->select_exist) {
				if ( (*iter)->subVar_exist ) 
					(*iter)->value = kNN_buffer[cnt].subject;
				else if ( (*iter)->geoVar_exist )
					(*iter)->value = kNN_buffer[cnt].geometry;
				else  
					(*iter)->value = kNN_buffer[cnt].restArgs[pos++]; 
			}
		}

   	observedOutputCardinality++;
   	cnt++;

   	return 1;
	}

#else

	#if RTREE_USED_IN_PLAN

		#if ONLY_POINT_GEOMS

			//cout << "next...WITHIN RTREE Selection Evaluation (POINT_GEOMS)" << endl;

			// Get the next tuple
			unsigned count = entry_RW->next();
			if (!count) return 0;

			observedOutputCardinality++;
			return 1;

		#else

			//cout << "next...WITHIN RTREE Selection Evaluation (ALL_GEOMS)" << endl;

			// Variables
			Result geoRes;
			double minX, maxX, minY, maxY;

			//cout << "HASGEO_ID: " << hasGeom_id_RW << endl;

			geoRes.flags = 0;  // reset the flags
			geoRes.flags = Result::idAvailable;

			// Get the next entry of priority queue 'SED'
			if (!SED.empty()) 
			{
				prRtreeEntry_RW = SED.top();
				SED.pop();

				#if LGD

					prRtreeEntry_RW.geometry = IS->getGeoID(prRtreeEntry_RW.subject, hasGeom_id_RW);  // get the ?geo id
					geoRes.id = prRtreeEntry_RW.geometry;

				#elif YAGO

					auto find_it = spatialIDs_TO_geomIDs.find(prRtreeEntry_RW.subject);

					if (find_it == spatialIDs_TO_geomIDs.end()) {  // not found
						prRtreeEntry_RW.geometry = IS->getGeoID(prRtreeEntry_RW.subject, hasGeom_id_RW);  // get the ?geo id
						geoRes.id = prRtreeEntry_RW.geometry;
						spatialIDs_TO_geomIDs.insert(make_pair(prRtreeEntry_RW.subject, prRtreeEntry_RW.geometry));	
					}
					else {  // found
						prRtreeEntry_RW.geometry = find_it->second;   
						geoRes.id = prRtreeEntry_RW.geometry;
					}

				#endif

				//cout << "GEOMETRY: " << prRtreeEntry_RW.geometry << endl;
			}
			else { return 0; }			

			geoRes.ensureString(this);  // get the ?geo string value ('retrieve geometry')

			type2_RW = getLongLat(geoRes.value, points2_RW, size2_RW);  // set 'points2_RW' & 'size2_RW' values (NORMALIZED COORDS)

   		// Get the MBR of the geometry
   		getRegion(points2_RW, size2_RW, minX, maxX, minY, maxY);

   		// Check the WITHIN predicate
   		if ( (minX >= points1_RW[0] && maxX <= points1_RW[2]) && (minY >= points1_RW[1] && maxY <= points1_RW[3]) )
   		{
      		//cout << "Pass Within WITH geometry retrieval" << endl;

      		// Pass register values to next operators	
      		unsigned pos = 0;
				for (vector<Register*>::const_iterator iter = entryTail_RW.begin(), limit = entryTail_RW.end(); iter != limit; ++iter) {
					if ((*iter)->select_exist) {
						if ( (*iter)->subVar_exist ) 
							(*iter)->value = prRtreeEntry_RW.subject;
						else if ( (*iter)->geoVar_exist )
							(*iter)->value = prRtreeEntry_RW.geometry;
						else  
							(*iter)->value = prRtreeEntry_RW.restArgs[pos++]; 
					}
				}

				observedOutputCardinality++;
				return 1;
   		}
   		else { 
   			cout << "False Positive is found" << endl;
   			return next(); 
   		}				

		#endif

	#else

   	//cout << "next...Selection Evaluation" << endl;

   	while (true) {  // till a result is produced

      	#if SEMIJOIN
      		unsigned count;
      		if (!done) {
	     	 		// Retrieve the next tuple
	      		count = input->next();
      		}
      		else count = 0;
      	#else
      		unsigned count = input->next();
	   		//cout << "id_left is "  << dynamic_cast<MergeJoin*>(input)->leftValue->value << endl;
	   		//cout << "id_right is " << dynamic_cast<MergeJoin*>(input)->rightValue->value << endl;
      	#endif

      	if (!count) 
      	{
				#if SEMIJOIN
					if (!done) {
						vIter = nonVerifiedEntries.begin();
						done = true;
						cout << "Verified entries: " << verifiedEntries.size() << endl;
						cout << "Non-verified entries: " << nonVerifiedEntries.size() << endl;
					}
					//check the remaining non-verified entries in the buffer
					for (;vIter!=nonVerifiedEntries.end();++vIter) {
		   			//cout << "vIter 1: " << vIter->first << " vIter 2: " << vIter->second << endl;
		   			string g1,g2;
		   			//retrieve geometries
		   			lookUpGeo(vIter->first,g1);
		   			lookUpGeo(vIter->second,g2);			   
		   			//std::cout << "Entered Distance selection..." << std::endl;
		   			//cout << "left coords: " << r1.value << " right coords: " << r2.value << endl;
		   			double n1,n2,n3,n4;
		   			//cout << "g1: " << g1 << " g2: " << g2 << endl;
		   			//get coordinates
		   			getLongLat(g1,n1,n2);
		   			getLongLat(g2,n3,n4);
		   			//cout << "n1: " << n1 << " n2: " << n2 << " n3: " << n3 << " n4: " << n4 << endl;
		   
		   			//XXX: hardcoded - only for <
		   			if (((n1-n3)*(n1-n3) + (n2-n4)*(n2-n4))<distThreshold) {
			 				//continue with the next entry
			 				unsigned checkedV = vIter->first;
			 				while (vIter->first==checkedV)
				 				++vIter;

			 				observedOutputCardinality+=1;
			 				return 1;
		   			}
					}
				#endif

				if (dynamic_cast<Selection::DistanceF*>(predicate)) {
					cout << "DistanceF. - FILTER" << endl;
					#if PRINTINFO
					std::cout << "- Total num of filtered tuples: " << filtered << endl;
					#endif
				}
				else if (dynamic_cast<Selection::WithinF*>(predicate)) {
					//cout << "WithinF. - FILTER" << endl;
					#if PRINTINFO
					std::cout << "- Total num of filtered tuples: " << filtered << endl;
					#endif
				}
				else if (dynamic_cast<Selection::Within*>(predicate)) {
					//cout << "Within. - NO FILTER" << endl;
					#if PRINTINFO
					std::cout << "Within Selection results: " << ((Selection::Within*) predicate)->cnt << std::endl;
					std::cout << "Within Selection double bindings: " << sum << std::endl;
					std::cout << "- Total num of verified tuples: " << verif << endl;
					#endif
				}
				else if (dynamic_cast<Selection::Distance*>(predicate)) {
					//cout << "Distance. - NO FILTER" << endl;
					//std::cout << "Selection results: " << ((Selection::Distance*)predicate)->cnt << std::endl;
					#if PRINTINFO
					std::cout << "Selection results: " << ((Selection::Distance*)predicate)->cnt << std::endl;
					std::cout << "Selection cache size: " << cache.size() << endl;
					std::cout << "Within Selection double bindings: " << sum << std::endl;
					std::cout << "- Total num of verified tuples: " << verif << endl;
					#endif
				}

				return 0;
      	}

      	// Match?
      	if (predicate->check()) {
      		//cout << "NEXT PASS" << endl;
	 			//if(count>1 && dynamic_cast<Selection::Distance*>(predicate)) std::cout << "count: " << count << endl;
        	 	observedOutputCardinality += count;
				//cout << "count: " << count << endl;
				sum += (count-1);
         	return count;
      	}
   	}

   #endif

#endif
}
//---------------------------------------------------------------------------
void Selection::print(PlanPrinter &out)
   // Print the operator tree. Debugging only
{
#if K_NEAREST_NEIGHBOR
	out.beginOperator("kNN_Selection", expectedOutputCardinality, observedOutputCardinality);
	entry->print(out);
	out.endOperator();
#else
	#if RTREE_USED_IN_PLAN
		out.beginOperator("Rtree_WITHIN_Selection", expectedOutputCardinality, observedOutputCardinality);
		entry_RW->print(out);
		out.endOperator();
	#else 
   	out.beginOperator("Selection", expectedOutputCardinality, observedOutputCardinality);
   	out.addGenericAnnotation(predicate->print(out));
   	input->print(out);
   	out.endOperator();
   #endif	
#endif
}
//---------------------------------------------------------------------------
void Selection::addMergeHint(Register *reg1, Register *reg2)
   // Add a merge join hint
{
#if K_NEAREST_NEIGHBOR
   entry->addMergeHint(reg1, reg2);
#else
   #if RTREE_USED_IN_PLAN
   	entry_RW->addMergeHint(reg1, reg2);
   #else	
		input->addMergeHint(reg1, reg2);
	#endif	
#endif
}
//---------------------------------------------------------------------------
void Selection::getAsyncInputCandidates(Scheduler &scheduler)
   // Register parts of the tree that can be executed asynchronous
{
#if K_NEAREST_NEIGHBOR
	entry->getAsyncInputCandidates(scheduler);
#else
	#if RTREE_USED_IN_PLAN
		entry_RW->getAsyncInputCandidates(scheduler);
	#else	
   	input->getAsyncInputCandidates(scheduler);
   #endif	
#endif
}
//---------------------------------------------------------------------------
void Selection::checkForGeo()
   // Register parts of the tree that can be executed asynchronous
{
#if K_NEAREST_NEIGHBOR
	entry->checkForGeo();
	isForGeo = entry->isForGeo;
#else
	#if RTREE_USED_IN_PLAN	
		entry_RW->checkForGeo();
		isForGeo = entry_RW->isForGeo;
	#else	
		input->checkForGeo();
		isForGeo = input->isForGeo;
	#endif	
#endif
}
//---------------------------------------------------------------------------