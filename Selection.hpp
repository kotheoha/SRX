#ifndef H_rts_operator_Selection
#define H_rts_operator_Selection
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
#include "rts/operator/Operator.hpp"
#include "rts/database/Database.hpp"
#include "rts/segment/FactsSegment.hpp"
#include "infra/util/Type.hpp"
#include <vector>
#include <string>
#include "cts/geo/Geometry.hpp"
#include <sstream>
#include "rts/cell/Grid.hpp"
#include "cts/infra/QueryGraph.hpp"
#include "cts/geo/GeoDef.hpp"
#include <unordered_map>
#include <unordered_set>
#include <queue>
using namespace std;
//---------------------------------------------------------------------------
class Register;
class Runtime;
//---------------------------------------------------------------------------
/// Applies a number of selections
class Selection : public Operator
{

   #if K_NEAREST_NEIGHBOR

      private:

      Operator *entry;  // the input
      vector<Register*> entryTail;  // all the attributes

      Geo::Geometry *geo;  // the geometry (it depends on encoding flag, if it IS or NOT in the Unit Square)
      Grid *grid;  // the grid

      uint32_t k_nearest_neighbor;  // the kNN k parameter value 
      bool kNN_encode;  // set 1 for Grid Encoding, 0 else
      double lastDist;    // the real distance of the 'k_st entity' from 'query'
      unsigned cell_cnt;  // number of processed 'cell entities' (extracting them from 'H')
      unsigned gridCells;  // total number of grid cells
      unsigned cnt, pos;
      
      double q_ll[2], q_ur[2], e_ll[2], e_ur[2];
      double rx1, ry1, rx2, ry2;
      double nx1, ny1, nx2, ny2;
      unsigned cX1, cY1, cX2, cY2;
      unsigned c_id, gpos, level;      
      double q_center[2], e_center[2];
      unsigned q_cell;

      double initDist[4];  // (U:0), (R:1), (D:2), (L:3) 

      double *points1, *points2;  // 'points1' is for the 'query geometry', 'points2' is for the 'retrieved geometries' 
      unsigned size1, size2;  // respective 'size' values of 'points1' and 'points2' 
      unsigned type1, type2;  // respective 'type' values of 'points1' and 'points2'

      struct kNNEntry {
         unsigned subject;
         unsigned geometry;
         bool bottomLevel; 
         unsigned zone;
         unsigned dir;  // (U:0), (R:1), (D:2), (L:3), (cell:4)
         double minDist;
         double realDist;
         vector<unsigned> restArgs;
      };

      kNNEntry priorityEntry;  // priority queue entry


      class CompareMinDistance {
         public:
         // d2 has higher priority than d1 if d2 priority field is smaller than d1
         bool operator() (kNNEntry &d1, kNNEntry &d2) {
            if (d2.minDist < d1.minDist) return true;
            else return false;   
         }
      };

      // MIN Distance priority queue
      priority_queue<kNNEntry, vector<kNNEntry>, CompareMinDistance> H;

      class CompareRealDistance {
         public:
         // d2 has higher priority than d1 if d2 priority field is bigger than d1
         bool operator() (kNNEntry &d1, kNNEntry &d2) {
            if (d2.realDist > d1.realDist) return true;
            else return false;   
         }
      };

      // REAL Distance priority queue
      priority_queue<kNNEntry, vector<kNNEntry>, CompareRealDistance> R;


      // BASELINE CASE < START > 
      struct kNN_res {
         unsigned subject;
         unsigned geometry;   
         double distance;
         vector<unsigned> restArgs;
      };

      vector<kNN_res> kNN_buffer;   

      struct sortByDistance {
         inline bool operator() (const kNN_res &lhs, const kNN_res &rhs) {
            return (lhs.distance < rhs.distance);
         }
      };
      // BASELINE CASE < END >

      // RTREE CASE FOR ALL GEOMS < START >  
      kNN_res prRtreeEntry;  // priority queue R-tree entry
      kNN_res prRtreeEntryTop;  // priority queue R-tree Top entry
   
      class CompareRealDistanceForRtree {
         public:
         // d2 has higher priority than d1 if d2 priority field is bigger than d1
         bool operator() (kNN_res &d1, kNN_res &d2) {
            if (d2.distance > d1.distance) return true;
            else return false;
         }
      };

      // REAL Distance priority queue
      priority_queue<kNN_res, vector<kNN_res>, CompareRealDistanceForRtree> RR;
      // RTREE CASE FOR ALL GEOMS < END >

   #endif


   private:

   #if YAGO   
      unordered_map<unsigned, unsigned> spatialIDs_TO_geomIDs;
   #endif

   struct WITHIN_res {
      unsigned subject;
      unsigned geometry;   
      vector<unsigned> restArgs;
   };

   WITHIN_res prRtreeEntry_RW;  // priority queue R-tree entry

   class CompareSpatialEntityIDsForRtree {
      public:
      // s2 has higher priority than s1 if s2 priority field is smaller than s1
      bool operator() (WITHIN_res &s1, WITHIN_res &s2) {
         if (s2.subject < s1.subject) return true;
         else return false;
      }
   };

   // Spatial Entity IDs priority queue
   priority_queue<WITHIN_res, vector<WITHIN_res>, CompareSpatialEntityIDsForRtree> SED;

   // RETRIEVE GEO IDs INDEX SCAN < START > 
   /// An index scan over the facts table
   class RetrieveGeoIDsIndexScan 
   {
      public:
      /// Hints during scanning
      class Hint : public FactsSegment::Scan::Hint {
         private:
         /// The scan
         RetrieveGeoIDsIndexScan& scan;

         public:
         /// Constructor
         Hint(RetrieveGeoIDsIndexScan& scan);
         /// Destructor
         ~Hint();
         /// The next hint
         void next(unsigned& value1,unsigned& value2,unsigned& value3);
      };
      friend class Hint;

      /// The facts segment
      FactsSegment& facts;
      /// The data order
      Database::DataOrder order;
      /// The scan
      FactsSegment::Scan scan;
      /// The hinting mechanism
      Hint hint;

      /// Constructor
      RetrieveGeoIDsIndexScan(Database& db, Database::DataOrder order);
      /// Destructor
      ~RetrieveGeoIDsIndexScan();
      /// Get the respective geo id based on filter1 & filter2 values 
      unsigned getGeoID(unsigned filter1, unsigned filter2);
   };
   // RETRIEVE GEO IDs INDEX SCAN < END >

   Operator *entry_RW;  // the input (RW: Rtree-WITHIN)
   vector<Register*> entryTail_RW;  // all the attributes
   Geo::Geometry *geo_RW;  // the geometry (in normalized coords)
   RetrieveGeoIDsIndexScan *IS;

   unsigned hasGeom_id_RW;  // the id of <hasGeometry> predicate

   double *points1_RW, *points2_RW;  // 'points1_RW' is for the 'query geometry', 'points2_RW' is for the 'retrieved geometries' 
   unsigned size1_RW, size2_RW;  // respective 'size' values of 'points1_RW' and 'points2_RW'
   unsigned type1_RW, type2_RW;  // respective 'type' values of 'points1_RW' and 'points2_RW'


   public:

   // a cache to avoid some dictionary lookups	
   unordered_map<unsigned, string> cache;
   unordered_map<unsigned, string>::iterator itter;

   /// A predicate result
   struct Result {
      /// Possible flags
      enum Flags { idAvailable = 1, stringAvailable = 2, typeAvailable = 4, subTypeAvailable = 8, booleanAvailable = 16 };

      /// The flags
      unsigned flags;
      /// The id
      unsigned id;
      /// The value
      std::string value;
      /// The type
      Type::ID type;
      /// The sub-type
      unsigned subType;
      /// The sub-type value
      std::string subTypeValue;
      /// The boolean interpretation
      bool boolean;

      /// Constructor
      Result() : flags(0) { }
      /// Destructor
      ~Result();

      /// Id available?
      bool hasId() const { return flags & idAvailable; }
      /// String available?
      bool hasString() const { return flags & stringAvailable; }

      /// Ensure that a string is available
      void ensureString(Selection *selection);
      /// Ensure that the type is available
      void ensureType(Selection *selection);
      /// Ensure that the subtype is available
      void ensureSubType(Selection *selection);
      /// Ensure that a boolean interpretation is available
      void ensureBoolean(Selection *selection);

      /// Set to a boolean value
      void setBoolean(bool v);
      /// Set to an id value
      void setId(unsigned v);
      /// Set to a string value
      void setLiteral(const std::string &c);
      /// Set to a string value
      void setIRI(const std::string &c);
   };

   /// Base for predicate evaluation
   class Predicate {
      protected:
      /// The outer selection
      Selection *selection;

      public:

      /// Constructor
      Predicate();
      /// Destructor
      virtual ~Predicate();

      /// Register the selection
      virtual void setSelection(Selection *selection);
      /// Evaluate the predicate
      virtual void eval(Result &result) = 0;
      /// Print the predicate (debugging only)
      virtual std::string print(PlanPrinter &out) = 0;

      /// Check the predicate
      bool check();
   };

   /// Binary operator
   class BinaryPredicate : public Predicate {
      protected:
      /// The input
      Predicate *left, *right;

      public:
      /// Constructor
      BinaryPredicate(Predicate *left, Predicate *right) : left(left), right(right) { }
      /// Destructor
      ~BinaryPredicate();

      /// Register the selection
      void setSelection(Selection *selection);
   };

   /// Unary operator
   class UnaryPredicate : public Predicate {
      protected:
      /// The input
      Predicate *input;

      public:
      /// Constructor
      UnaryPredicate(Predicate *input) : input(input) { }
      /// Destructor
      ~UnaryPredicate();

      /// Register the selection
      void setSelection(Selection *selection);
   };

   /// WithinF selection
   class WithinF : public UnaryPredicate {
      double c1[2], c2[2], rx1, ry1, rx2, ry2;
      bool res;
      unsigned cX1, cY1, cX2, cY2;
      unsigned c_id, gpos, level;

      public:
      Grid *grid;
      Geo::Geometry *geo;  // the normalized geometry (IN the Unit Square)

      /// Constructor
      WithinF(Predicate *input, Geo::Geometry *g, Grid *gr) : UnaryPredicate(input), geo(new Geo::Rectangle((Geo::Rectangle*) g)), grid(gr) {
			geo->toUnitSquare(grid->minX, grid->minY, grid->maxX, grid->maxY);
			//cout << "Geometry for WithinF selection: " << endl;
			//geo->print();
			//cout << endl;
		}

      /// Evaluate the predicate
      void eval(Result &result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter &out);

		// get original points
		void getUpperRightPoint(unsigned x, unsigned y, double *coord);
		void getLowerLeftPoint(unsigned x, unsigned y, double *coord);

		// get c_id, g_pos and level
		bool getCellId(unsigned value, unsigned &c_id, unsigned &g_pos, unsigned &level); 
   };

   /// DistanceF selection
   class DistanceF : public Predicate {
      private:
	   /// Arguments (cell ids)
	   Predicate *arg1, *arg2, *arg3;

		bool res;
		Grid *grid;
		double SQRT_2;

      // the square of the distance (NOT in the Unit Square) and its x,y offsets in the Unit Square
      double distThreshold, minDist, maxDist;	

		// the type of join
		QueryGraph::Filter::Type type;

		// get cell ids
		unsigned lcID, level1, g_pos1, rcID, level2, g_pos2;
		unsigned leftCellX1, leftCellY1, leftCellX2, leftCellY2, rightCellX1, rightCellY1, rightCellX2, rightCellY2;

		// get c_id, g_pos and level
		bool getCellId(unsigned value, unsigned &c_id, unsigned &g_pos, unsigned &level);

		// get original points
		void getUpperRightPoint(unsigned x, unsigned y, double *coord);
		void getLowerRightPoint(unsigned x, unsigned y, double *coord);
		void getUpperLeftPoint(unsigned x, unsigned y, double *coord);
		void getLowerLeftPoint(unsigned x, unsigned y, double *coord);

		// compute distance^2
		double dist2(double *c1, double *c2);

      void getCenterCoords(unsigned cord_x1, unsigned cord_y1, unsigned cord_x2, unsigned cord_y2, double *c);
	   bool evalPredicate(double minDist, double maxDist, bool &verified);
   	void getMinMaxDist(unsigned cord_x1, unsigned cord_y1, unsigned cord_x2, unsigned cord_y2, unsigned c_x1, 
								 unsigned c_y1, unsigned c_x2, unsigned c_y2, double &minDist, double &maxDist);
		double param1, param2;

      public:

      /// Constructor
      DistanceF(Predicate *arg1, Predicate *arg2, Predicate *arg3, QueryGraph::Filter::Type jType, Grid *g) : 
		arg1(arg1), arg2(arg2), arg3(arg3), grid(g), type(jType)
      {
			// cast the distance threshold only once in the beginning
      	stringstream ss(stringstream::in | stringstream::out);
			Result r;
			arg3->eval(r);
			ss << r.value;
			ss >> distThreshold;
			distThreshold *= distThreshold;  // e^2
         //cout << "DistanceF threshold: " << distThreshold << endl;

			SQRT_2 = sqrt(2);

			param1 = (grid->cSide * grid->lengthX) * (grid->cSide * grid->lengthX);
			param2 = (grid->cSide * grid->lengthY) * (grid->cSide * grid->lengthY);
      };
      /// Destructor
      ~DistanceF();

      /// Register the selection
      void setSelection(Selection *selection);
      /// Evaluate the predicate
      void eval(Result &result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter &out);
   };

   #ifdef HASGEO
   /// Within selection
   class Within : public UnaryPredicate {
      public:

      Geo::Geometry *geo, *geo2;  // the normalized geometry (NOT in the Unit Square)
      double minX, maxX, minY, maxY;
      double *points;
      unsigned cnt, size;
    
      Grid *gr;

      /// Constructor
      Within(Predicate *input, Geo::Geometry *g, Grid *grid) : UnaryPredicate(input), geo(new Geo::Rectangle((Geo::Rectangle*) g)), gr(grid)
	   {
		   #if GEOSTORE
		      geo2 = new Geo::Rectangle((Geo::Rectangle*) geo);
			   geo2->toUnitSquare(gr->minX, gr->minY, gr->maxX, gr->maxY);
			   geo->normalize();	
	      #else
			   geo->normalize();
		   #endif
		   cnt = 0;
	   }

      /// Evaluate the predicate
      void eval(Result &result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter &out);
      void getLongLat(const std::string &geometry, double* &points, unsigned &size);
      void getRegion(double *points, unsigned size, double &minX, double &maxX, double &minY, double &maxY);
      unsigned countPoints(const string &geom, const string &dlm);
   };

   /// Distance selection
   class Distance : public BinaryPredicate {
      public:
		double distThreshold;	//the square of the distance 
		double n1,n2,n3,n4;
		double *points1,*points2;

		//the type of join
		QueryGraph::Filter::Type type;
      unsigned cnt, size1,size2;

      /// Constructor
      Distance(Predicate* left,Predicate* right, QueryGraph::Filter::Type jType,const std::string& dist2) : BinaryPredicate(left,right),type(jType)
		{
			//cast the distance threshold only once in the beginning
			stringstream ss(stringstream::in|stringstream::out);
			ss << dist2;
			ss >> distThreshold;
	
			distThreshold *= distThreshold;		//e^2

			//cout << "Distance threshold: " << distThreshold << endl;

			cnt = 0;
		}
      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);

      unsigned countPoints(const string& geom, const string& dlm);
      unsigned getLongLat(const std::string& geometry, double*& points, unsigned& size);
      void getRegion(double *points, unsigned size, double& minX, double& maxX, double& minY, double& maxY);

      double dot(double* A, double* B);
      double distance2(double* A, double* B);
      double pointLineDist2(double* point, double* segment);
      double lineLineDist2(double* seg1, double* seg2);
      double linePolygonDist2(double* line, double* polygon, unsigned size);
      double pointPolygonDist2(double* point, double* polygon, unsigned size);
      double polygonPolygonDist2(double* polygon1, double* polygon2, unsigned size1, unsigned size2);
      double pointMultipointDist2(double* point, double* multipoint, unsigned size);
      double lineMultipointDist2(double* line, double* multipoint, unsigned size);
      double polygonMultipointDist2(double* polygon, double* multipoint, unsigned size1, unsigned size2);
      double multipointMultipointDist2(double* multipoint1, double* multipoint2, unsigned size1, unsigned size2);

      double pointLineStringDist2(double* point, double* line, unsigned size);
      double lineStringLineStringDist2(double* line1, unsigned size1, double* line2, unsigned size2);
      double lineStringPolygonDist2(double* line, unsigned size1, double* polygon, unsigned size2);
      double lineStringMultipointDist2(double* line, unsigned size1, double* multipoint, unsigned size2);
   };
   #endif

   /// Logical or
   class Or : public BinaryPredicate {
      public:
      /// Constructor
      Or(Predicate* left,Predicate* right) : BinaryPredicate(left,right) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Logical and
   class And : public BinaryPredicate {
      public:
      /// Constructor
      And(Predicate* left,Predicate* right) : BinaryPredicate(left,right) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Comparison ==
   class Equal : public BinaryPredicate {
      public:
      /// Constructor
      Equal(Predicate* left,Predicate* right) : BinaryPredicate(left,right) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Comparison !=
   class NotEqual : public BinaryPredicate {
      public:
      /// Constructor
      NotEqual(Predicate* left,Predicate* right) : BinaryPredicate(left,right) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Comparison <
   class Less : public BinaryPredicate {
      public:
      /// Constructor
      Less(Predicate* left,Predicate* right) : BinaryPredicate(left,right) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Comparison <=
   class LessOrEqual : public BinaryPredicate {
      public:
      /// Constructor
      LessOrEqual(Predicate* left,Predicate* right) : BinaryPredicate(left,right) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Arithmetic +
   class Plus : public BinaryPredicate {
      public:
      /// Constructor
      Plus(Predicate* left,Predicate* right) : BinaryPredicate(left,right) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Arithmetic -
   class Minus : public BinaryPredicate {
      public:
      /// Constructor
      Minus(Predicate* left,Predicate* right) : BinaryPredicate(left,right) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Arithmetic *
   class Mul : public BinaryPredicate {
      public:
      /// Constructor
      Mul(Predicate* left,Predicate* right) : BinaryPredicate(left,right) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Arithmetic /
   class Div : public BinaryPredicate {
      public:
      /// Constructor
      Div(Predicate* left,Predicate* right) : BinaryPredicate(left,right) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Operator !
   class Not : public UnaryPredicate {
      public:
      /// Constructor
      Not(Predicate* input) : UnaryPredicate(input) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Operator -
   class Neg : public UnaryPredicate {
      public:
      /// Constructor
      Neg(Predicate* input) : UnaryPredicate(input) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// A NULL value
   class Null : public Predicate {
      public:
      /// Evaluate the predicate
      void eval(Result &result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter &out);
   };

   /// A false value
   class False : public Predicate {
      public:
      /// Evaluate the predicate
      void eval(Result &result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter &out);
   };

   /// Variable access
   class Variable : public Predicate {
      protected:
      /// The register
      Register *reg;

      public:
      /// Constructor
      Variable(Register *reg) : reg(reg) { }

      void setResolved(bool f);
      /// Evaluate the predicate
      void eval(Result &result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter &out);
   };

   /// Constant
   class ConstantLiteral : public Predicate {
      private:
      /// The id
      unsigned id;

      public:
      /// Constructor
      ConstantLiteral(unsigned id) : id(id) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Constant
   class TemporaryConstantLiteral : public Predicate {
      private:
      /// The value
      std::string value;

      public:
      /// Constructor
      TemporaryConstantLiteral(const std::string& value) : value(value) {}

      /// Evaluate the predicate
      void eval(Result& result);
      //return the actual value in the form of string
      std::string getValue();
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Constant
   class ConstantIRI : public Predicate {
      private:
      /// The id
      unsigned id;

      public:
      /// Constructor
      ConstantIRI(unsigned id) : id(id) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Constant
   class TemporaryConstantIRI : public Predicate {
      private:
      /// The value
      std::string value;

      public:
      /// Constructor
      TemporaryConstantIRI(const std::string& value) : value(value) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Function call
   class FunctionCall : public Predicate {
      private:
      /// The function
      std::string func;
      /// Arguments
      std::vector<Predicate*> args;

      public:
      /// Constructor
      FunctionCall(const std::string& func,const std::vector<Predicate*>& args) : func(func),args(args) {}

      /// Register the selection
      void setSelection(Selection* selection);

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Builtin str
   class BuiltinStr : public UnaryPredicate {
      public:
      /// Constructor
      BuiltinStr(Predicate* input) : UnaryPredicate(input) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Builtin lang
   class BuiltinLang : public UnaryPredicate {
      public:
      /// Constructor
      BuiltinLang(Predicate* input) : UnaryPredicate(input) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Builtin langMatches
   class BuiltinLangMatches : public BinaryPredicate {
      public:
      /// Constructor
      BuiltinLangMatches(Predicate* left,Predicate* right) : BinaryPredicate(left,right) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Builtin datatype
   class BuiltinDatatype : public UnaryPredicate {
      public:
      /// Constructor
      BuiltinDatatype(Predicate* input) : UnaryPredicate(input) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Builtin bound
   class BuiltinBound : public Predicate {
      private:
      /// The register
      Register* reg;

      public:
      /// Constructor
      BuiltinBound(Register* reg) : reg(reg) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Builtin sameTerm
   class BuiltinSameTerm : public BinaryPredicate {
      public:
      /// Constructor
      BuiltinSameTerm(Predicate* left,Predicate* right) : BinaryPredicate(left,right) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Builtin isIRI
   class BuiltinIsIRI : public UnaryPredicate {
      public:
      /// Constructor
      BuiltinIsIRI(Predicate* input) : UnaryPredicate(input) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Builtin isBlank
   class BuiltinIsBlank : public UnaryPredicate {
      public:
      /// Constructor
      BuiltinIsBlank(Predicate* input) : UnaryPredicate(input) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Builtin isLiteral
   class BuiltinIsLiteral : public UnaryPredicate {
      public:
      /// Constructor
      BuiltinIsLiteral(Predicate* input) : UnaryPredicate(input) {}

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Builtin RegEx
   class BuiltinRegEx : public Predicate {
      private:
      /// Arguments
      Predicate* arg1,*arg2,*arg3;

      public:
      /// Constructor
      BuiltinRegEx(Predicate* arg1,Predicate* arg2,Predicate* arg3) : arg1(arg1),arg2(arg2),arg3(arg3) {}
      /// Destructor
      ~BuiltinRegEx();

      /// Register the selection
      void setSelection(Selection* selection);

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   /// Builtin in
   class BuiltinIn : public Predicate {
      private:
      /// The probe
      Predicate* probe;
      /// The set
      std::vector<Predicate*> args;

      public:
      /// Constructor
      BuiltinIn(Predicate* probe,const std::vector<Predicate*>& args) : probe(probe),args(args) {}

      /// Register the selection
      void setSelection(Selection* selection);

      /// Evaluate the predicate
      void eval(Result& result);
      /// Print the predicate (debugging only)
      std::string print(PlanPrinter& out);
   };

   protected:
   /// The input
   Operator *input;
   /// The runtime
   Runtime &runtime;
   /// The predicate
   Predicate *predicate;

   unsigned sum;

	#if SEMIJOIN
	   unordered_set<unsigned> verifiedEntries;
		multimap<unsigned,unsigned> nonVerifiedEntries;
		multimap<unsigned,unsigned>::iterator vIter,vIterStop;
	
		void lookUpGeo(unsigned id, string& g);
		void getLongLat(string& g, double& lon, double& lat);

		bool done;

		double distThreshold;
	#endif

   public:

   // Get the coordinates of the MBR for a given geometry   
   void getRegion(double *points, unsigned size, double &minX, double &maxX, double &minY, double &maxY);
   // Return the type of geometry (0:point, 1:line, 2:polygon, 3:multipoint)
   unsigned getLongLat(const string &geometry, double* &points, unsigned &size);

   /// Constructor
   Selection(Operator *input, Runtime &runtime, Predicate *predicate, double expectedOutputCardinality);

   /// Constructor to check false positives when R-tree is used to evaluate the WITHIN predicate
   Selection(Operator *entry_RW, vector<Register*> &entryTail_RW, Runtime &runtime,
             Geo::Geometry *g, double expectedOutputCardinality);

	#if K_NEAREST_NEIGHBOR
		/// Constructor kNN
		Selection(Operator *entry, vector<Register*> &entryTail, Runtime &runtime, Geo::Geometry *g, 
				 	 uint32_t k, double expectedOutputCardinality, Grid *gr, bool encodeFlag); 

		unsigned long long getHilbertId(double normX1, double normX2, double normY1, double normY2, double normCellSide, 
												  unsigned **c2i, unsigned cellsPerDim, unsigned b1, unsigned b2, unsigned &level);
		unsigned getGridPos(unsigned cId, unsigned level, unsigned size, unsigned MAX_LEVEL, unsigned b1);
		bool getCellId(unsigned value, unsigned &c_id, unsigned &g_pos, unsigned &level);

      // Get the geometry and assign points and size values
		void parseGeometry(Geo::Geometry *g, double* &points, unsigned &size);

      // Get the minimum distance between a point and a region based on libspatialindex library definition
      double RtreeDistPointRegion(double *point, double *region, unsigned size);

		// Find the 'limit entity conditions' for each level based on 'maxCell'
		void findEntityLimits(unsigned *limits, unsigned* &grid, unsigned maxCell, unsigned b2);
		// Find the maximum bottom cell id of '(zone, dir)' rectangle zone	
		unsigned findMaxCellInRecZone(unsigned zone, unsigned dir, unsigned** &c2i, double *q_ll, double cSide);

      double dot(double *A, double *B);
      double distance2(double *A, double *B);
      double pointLineDist2(double *point, double *segment);
		double pointLineStringDist2(double *point, double *line, unsigned size);
      double lineLineDist2(double *seg1, double *seg2);
		double lineStringLineStringDist2(double *line1, unsigned size1, double *line2, unsigned size2);
      double linePolygonDist2(double *line, double *polygon, unsigned size);
		double lineStringPolygonDist2(double *line, unsigned size1, double *polygon, unsigned size2);
      double pointPolygonDist2(double *point, double *polygon, unsigned size);
      double polygonPolygonDist2(double *polygon1, double *polygon2, unsigned size1, unsigned size2);
      double pointMultipointDist2(double *point, double *multipoint, unsigned size);
      double lineMultipointDist2(double *line, double *multipoint, unsigned size);
		double lineStringMultipointDist2(double *line, unsigned size1, double *multipoint, unsigned size2);	
      double polygonMultipointDist2(double *polygon, double *multipoint, unsigned size1, unsigned size2);
      double multipointMultipointDist2(double *multipoint1, double *multipoint2, unsigned size1, unsigned size2); 
	#endif

   /// Destructor
   ~Selection();

   /// Produce the first tuple
   unsigned first();
   /// Produce the next tuple
   unsigned next();

   void checkForGeo();

   /// Print the operator tree. Debugging only.
   void print(PlanPrinter &out);
   /// Add a merge join hint
   void addMergeHint(Register *reg1, Register *reg2);
   /// Register parts of the tree that can be executed asynchronous
   void getAsyncInputCandidates(Scheduler &scheduler);
};
//---------------------------------------------------------------------------
#endif
