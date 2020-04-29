#ifndef H_rts_runtime_BulkOperation
#define H_rts_runtime_BulkOperation
//---------------------------------------------------------------------------
#include "rts/runtime/DifferentialIndex.hpp"
#include "rts/runtime/PredicateLockManager.hpp"
#include <map>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
//---------------------------------------------------------------------------
// RDF-3X
// (c) 2009 Thomas Neumann. Web site: http://www.mpi-inf.mpg.de/~neumann/rdf3x
//
// This work is licensed under the Creative Commons
// Attribution-Noncommercial-Share Alike 3.0 Unported License. To view a copy
// of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/
// or send a letter to Creative Commons, 171 Second Street, Suite 300,
// San Francisco, California, 94105, USA.
//---------------------------------------------------------------------------
/// A bulk operation/transaction
class BulkOperation
{
   private:

   /// The differential index
   DifferentialIndex &differentialIndex;
	/// The triples
	std::vector<DifferentialIndex::Triple> triples;	

	/// The triples for 'delete' from the database indexes	
	std::set<DifferentialIndex::SetTriple> triples_for_delete;
	/// The triples for 'insert' to the database indexes
	std::set<DifferentialIndex::SetTriple> triples_for_insert;
	/// The triples for 'insert' to the database indexes (2)
	std::set<DifferentialIndex::SetTriple> triples_for_insert_2;

   /// The temporary dictionary
   std::map<DifferentialIndex::Literal, unsigned> string2id;
   /// The temporary dictionary
   std::vector<DifferentialIndex::Literal> id2string;
   /// Will we delete entries?
   bool deleteMarker;

	// Its entry can be either 
	// 'the new ID of an entity that changed its ID due to updates input'
 	// or
	// 'the ID of a non-exist(new) to database entity that extends the dictionary' 
	std::unordered_set<unsigned> notInDB;

	// Used for Batch Statistics < START >
	// Batch ID
	unsigned cntB;	
	// Database Lookups Counter in Filter_Stage
	unsigned long fL_DB_C;
	// Dictionary Lookups Counter in Filter_Stage
	unsigned long fL_dict_C;
	// Number of Triples Inserted in Facts B-Trees	
	unsigned long tI_C;	
	// Number of Triples Deleted From Facts B-Trees
	unsigned long tD_C;
	// Number of ReEncodings
	unsigned long ReEncode_C;
	//unsigned long tI_ID_C, tI_DEL_C;
	// Used for Batch Statistics < END >

   /// Map a string
   unsigned mapString(const DifferentialIndex::Literal &value);

   public:

   /// Constructor
   explicit BulkOperation(DifferentialIndex &differentialIndex);
   /// Destructor
   ~BulkOperation();

   /// Add a triple
   void insert(const std::string &subject, const std::string &predicate, const std::string &object,
					Type::ID objectType, const std::string &objectSubType);
	/// Add a triple (updates version)
   void insert_updates(const std::string &subject, const std::string &predicate, const std::string &object,
					        Type::ID objectType, const std::string &objectSubType, bool op);
   /// Build a lock cover
   void buildCover(unsigned maxSize, std::vector<PredicateLockManager::Box> &boxes);
   /// Delete operation
   void markDeleted() { deleteMarker = true; }
	/// Initialize the parameters of a batch
	void initBatchParams(unsigned bID) { 
		cntB = bID;  fL_DB_C = 0;  fL_dict_C = 0;  tI_C = 0;  tD_C = 0;  ReEncode_C = 0;
	}	
   /// Commit
   void commit();
	/// Commit (updates version)
	void commit_updates(std::fstream &f2);
   /// Abort
   void abort();	
};
//---------------------------------------------------------------------------
#endif
