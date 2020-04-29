#include "rts/runtime/BulkOperation.hpp"
#include "rts/database/Database.hpp"
#include "rts/segment/DictionarySegment.hpp"
#include "infra/util/Type.hpp"
#include <algorithm>
#include <cassert>
#include "rts/segment/FactsSegment.hpp"
#include "rts/segment/FullyAggregatedFactsSegment.hpp"
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
using namespace std;
//---------------------------------------------------------------------------
BulkOperation::BulkOperation(DifferentialIndex &differentialIndex)
   : differentialIndex(differentialIndex), deleteMarker(false)
   // Constructor
{
}
//---------------------------------------------------------------------------
BulkOperation::~BulkOperation()
   // Destructor
{
}
//---------------------------------------------------------------------------
unsigned BulkOperation::mapString(const DifferentialIndex::Literal& value)
   // Map a string
{
   // Local?
   if (string2id.count(value)) 
      return string2id[value];

   // Resolve the sub-type if any
   unsigned subType=0;
   if (Type::hasSubType(value.type)) {
      if (value.type==Type::CustomType) {
         Type::ID realType=value.type;
         if (value.subType=="http://www.w3.org/2001/XMLSchema#string") {
            realType=Type::String;
         } else if (value.subType=="http://www.w3.org/2001/XMLSchema#integer") {
            realType=Type::Integer;
         } else if (value.subType=="http://www.w3.org/2001/XMLSchema#decimal") {
            realType=Type::Decimal;
         } else if (value.subType=="http://www.w3.org/2001/XMLSchema#double") {
            realType=Type::Double;
         } else if (value.subType=="http://www.w3.org/2001/XMLSchema#boolean") {
            realType=Type::Boolean;
         }
         if (realType!=value.type) {
            DifferentialIndex::Literal l;
            l.value=value.value;
            l.type=realType;
            return mapString(l);
         }
      }
      DifferentialIndex::Literal l;
      l.value=value.subType;
      l.type=Type::getSubTypeType(value.type);
      subType=mapString(l);
   };

   // Already in db?
   unsigned id;
   if (differentialIndex.lookup(value.value,value.type,subType,id))
      return id;

   // Create a temporary id
   id=(~0u)-id2string.size();
   string2id[value]=id;
   id2string.push_back(value);

   return id;
}
//---------------------------------------------------------------------------
void BulkOperation::insert(const string &subject, const string &predicate, const string &object,
								   Type::ID objectType, const std::string &objectSubType)
	// Add a triple
{
	DifferentialIndex::Triple t;
	DifferentialIndex::Literal l;	

   l.value = subject;
   l.type = Type::URI;
   t.subject = mapString(l);

   l.value = predicate;
   t.predicate = mapString(l);

   l.value = object;
   l.type = objectType;
   l.subType = objectSubType;
   t.object = mapString(l);

	triples.push_back(t);
}
//---------------------------------------------------------------------------
void BulkOperation::insert_updates(const string &subject, const string &predicate, const string &object,
								           Type::ID objectType, const string &objectSubType, bool op)
   // Add a triple (updates version)
{
	#if UPDATE_LOG
	cout << "'BO::insert(..)' FUNCTION" << endl;
	cout << "s_str: " << subject << ", p_str: " << predicate << ", o_str: " << object 
		  << ", o_str_type: " << objectType << ", o_str_subType: " << objectSubType << " op: " << op << endl; 
	#endif	
	

#if 0
	DifferentialIndex::SetTriple t;	
	string hasGeo = "hasGeometry";

	Type::ID objectSubType_type;
	unsigned objectSubType_id;

	// get the subject ID 't.subject'
	differentialIndex.getDictionary().string2id_updates(subject, Type::URI, 0, t.subject);
	// get the predicate ID 't.predicate'
	differentialIndex.getDictionary().string2id_updates(predicate, Type::URI, 0, t.predicate);
	// get the object ID 't.object'
	if (Type::hasSubType(objectType)) {	
		objectSubType_type = Type::getSubTypeType(objectType);
		differentialIndex.getDictionary().string2id_updates(objectSubType, objectSubType_type, 0, objectSubType_id);
		differentialIndex.getDictionary().string2id_updates(object, objectType, objectSubType_id, t.object);
	}
	else differentialIndex.getDictionary().string2id_updates(object, Type::URI, 0, t.object); 

	char *eValue;
	unsigned eValueSize;
	Type::ID type;
	unsigned subType; 

	differentialIndex.getDictionary().id2string_updates(t.subject, eValue, eValueSize, type, subType);
	string SUB = eValue; 
	cout << "SUB: " << SUB << endl; cout << "SUB_LEN_1: " << SUB.length() << endl; cout << "SUB_LEN_2: " << eValueSize << endl;
	cout << "SUB_TYPE: " << type << endl; cout << "SUB_SUBTYPE: " << subType << endl;   
		
	delete[] eValue; cout << endl;

	differentialIndex.getDictionary().id2string_updates(t.predicate, eValue, eValueSize, type, subType);
	string PRED = eValue;
	cout << "PRED: " << PRED << endl; cout << "PRED_LEN_1: " << PRED.length() << endl; cout << "PRED_LEN_2: " << eValueSize << endl;
	cout << "PRED_TYPE: " << type << endl; cout << "PRED_SUBTYPE: " << subType << endl; 

	delete[] eValue; cout << endl;

	differentialIndex.getDictionary().id2string_updates(t.object, eValue, eValueSize, type, subType);
	if (predicate.compare(hasGeo) == 0) {
		unsigned geoType = *(unsigned*) eValue;
		double *points = (double*) &eValue[sizeof(unsigned)]; 
		unsigned size = (eValueSize - sizeof(unsigned))/sizeof(double);
		cout << "GEO_TYPE: " << geoType << endl; cout << "size: " << size << endl;
		for (unsigned i = 0; i < size; i = i + 2)
			cout << "(" << points[i] << ", " << points[i+1] << ")" << " ";
		cout << endl;	
	}  
	else {
		string OBJ = eValue;
		cout << "OBJ: " << OBJ << endl; cout << "OBJ_LEN_1: " << OBJ.length() << endl; cout << "OBJ_LEN_2: " << eValueSize << endl;
		cout << "OBJ_TYPE: " << type << endl; cout << "OBJ_SUBTYPE: " << subType << endl; 
	}

	delete[] eValue; cout << endl;
#endif


#if 1
	DifferentialIndex::SetTriple t;
	string hasGeo = "hasGeometry";

	Type::ID objectSubType_type;
	unsigned objectSubType_id;

	FactsSegment::Scan scan_SPO, scan_OPS;

	// We assume that all 'triple SUBJECTS' have type URI
	differentialIndex.getDictionary().string2id_updates(subject, Type::URI, 0, t.subject);  fL_dict_C++;
	
	// CASE: 's' NOT in Dictionary < START >
	if (t.subject == ~0u) 
	{	
		#if UPDATE_LOG
		cout << "** UPDATES (1/3) CASE **" << endl;
		#endif
		// FIRST WAY...	
		if ( (op == true) && (predicate.compare(hasGeo) == 0) ) 
		{		
			#if UPDATE_LOG	
			cout << "FIRST WAY" << endl;
			#endif	
			differentialIndex.getDictionary().updateDictWithGrid(
				subject, predicate, object, objectType, objectSubType, differentialIndex.getDatabase().grid, 
				1, t.subject, t.predicate, t.object, fL_dict_C
			);	
			// skip duplicate <hasGeometry> triples
			if (!( (t.subject == 0) && (t.predicate == 0) && (t.object == 0) )) {
				triples_for_insert.insert(t);
				notInDB.insert(t.subject);	
			}
		}
		// SECOND WAY...
		else if ( (op == true) && (predicate.compare(hasGeo) != 0) ) 
		{
			#if UPDATE_LOG
			cout << "SECOND WAY" << endl;
			#endif
			differentialIndex.getDictionary().updateDictionary(t.subject, 0, subject, Type::URI, 0);  fL_dict_C++;
			// We assume that all 'triple PREDICATES' have type URI
			differentialIndex.getDictionary().string2id_updates(predicate, Type::URI, 0, t.predicate);  fL_dict_C++;
			if (t.predicate == ~0u) {
				#if UPDATE_LOG	
				cout << "PRED: " << predicate << " NOT in Dict" << endl;
				#endif
				differentialIndex.getDictionary().updateDictionary(t.predicate, 0, predicate, Type::URI, 0);  fL_dict_C++;
				notInDB.insert(t.predicate);	
			} 
			// We assume that all 'triple OBJECTS' have two separate cases
			if (Type::hasSubType(objectType)) {	
				#if UPDATE_LOG
				cout << "OBJ: " << object << " is NOT a simple string (Type != URI)" << endl;
				#endif
				objectSubType_type = Type::getSubTypeType(objectType);
				differentialIndex.getDictionary().string2id_updates(objectSubType, objectSubType_type, 0, objectSubType_id);  fL_dict_C++;
				if (objectSubType_id == ~0u) {
					#if UPDATE_LOG
					cout << "OBJ_SubType: " << objectSubType << " NOT in Dict" << endl;	
					#endif
					differentialIndex.getDictionary().updateDictionary(objectSubType_id, 0, objectSubType, objectSubType_type, 0);  fL_dict_C++;	
					differentialIndex.getDictionary().updateDictionary(t.object, 0, object, objectType, objectSubType_id);  fL_dict_C++;			
					notInDB.insert(t.object);		
				}	
				else {
					differentialIndex.getDictionary().string2id_updates(object, objectType, objectSubType_id, t.object);  fL_dict_C++;
					if (t.object == ~0u) {	
						#if UPDATE_LOG			
						cout << "OBJ: " << object << " NOT in Dict" << endl;
						#endif
						differentialIndex.getDictionary().updateDictionary(t.object, 0, object, objectType, objectSubType_id);  fL_dict_C++;
						notInDB.insert(t.object);
					}
				}
			}
			else { 
				#if UPDATE_LOG	
				cout << "OBJ: " << object << " is a simple string (Type = URI or Literal)" << endl;
				#endif
				differentialIndex.getDictionary().string2id_updates(object, objectType, 0, t.object);  fL_dict_C++;
				if (t.object == ~0u) {
					#if UPDATE_LOG						
					cout << "OBJ: " << object << " NOT in Dict" << endl;
					#endif
					differentialIndex.getDictionary().updateDictionary(t.object, 0, object, objectType, 0);  fL_dict_C++;
					notInDB.insert(t.object);
				}
			}
			triples_for_insert.insert(t);
			triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
			notInDB.insert(t.subject);
		}
	}
	// CASE: 's' NOT in Dictionary < END >

 	// ............................................................................... //

	// CASE: 's' EXISTS in Dictionary AND 's' is NON-SPATIAL < START > 
	else if ( (t.subject != ~0u) && ((t.subject % 2) == 0) ) 
	{
		#if UPDATE_LOG	
		cout << "** UPDATES (2/3) CASE **" << endl;
		#endif
		// FIRST WAY...
		if ( (op == true) && (predicate.compare(hasGeo) == 0) ) 
		{
			#if UPDATE_LOG
			cout << "FIRST WAY..." << endl;
			#endif	
			unsigned s_old = t.subject;  // the subject value that will change
			differentialIndex.getDictionary().updateDictWithGrid(
				subject, predicate, object, objectType, objectSubType, differentialIndex.getDatabase().grid,
				0, t.subject, t.predicate, t.object, fL_dict_C
			);			
			// skip duplicate <hasGeometry> triples
			if (!( (t.subject == 0) && (t.predicate == 0) && (t.object == 0) )) 
			{
				unsigned s_new = t.subject;  // the changed ID

				triples_for_insert.insert(t);	

				#if UPDATE_LOG
				cout << "(SUB_OLD, SUB_NEW): (" << s_old << ", " << s_new << ")" << endl; 
				#endif
				differentialIndex.getDictionary().updateDictForIDs(subject, Type::URI, 0, s_old, s_new);  fL_dict_C++; 
				ReEncode_C++;
				
				// CASE_1: "s_old" NOT in Database
				if (notInDB.count(s_old))
				{
					#if UPDATE_LOG
					cout << "CASE_1: SIMPLE" << endl;
					#endif	

					// Update all the 's_old-triples' with the new s value 's_new' in   	
					// both 'triples_for_insert' and 'triples_for_insert_2' triple sets  
					// (we consider safe to not examine the 'predicate' triple places) < START >	 
					// ..........................................................................
					vector<DifferentialIndex::SetTriple> forMain, forSecondary;
					vector<DifferentialIndex::SetTriple> forInsert;

					// "Note": for (ub1,ub2) and (ub3,ub4) it holds that if 's_old' does not exist 
					// in 'triples_for_insert' and 'triples_for_insert_2' as subject within a triple, 
					// then ub1 = ub2 and ub3 = ub4, so the respective 'for' will not be executed  

					// points to the first triple that starts with 's_old' as subject in 'triples_for_insert'
					auto ub1 = triples_for_insert.upper_bound(DifferentialIndex::SetTriple(s_old,0,0));
					// points to the first triple that starts with '> s_old' as subject in 'triples_for_insert'
					auto ub2 = triples_for_insert.upper_bound(DifferentialIndex::SetTriple((s_old+1),0,0));
					for (auto it = ub1; it != ub2; it++) {
						unsigned v2, v3;
						v2 = it->predicate; v3 = it->object;
						//#if (!LGD)
						forSecondary.push_back(DifferentialIndex::SetTriple(v3,v2,s_old));
						//#endif	
						forInsert.push_back(DifferentialIndex::SetTriple(s_new,v2,v3));
					}		
					// Replacement
					if (ub1 != ub2) triples_for_insert.erase(ub1,ub2); 
					for (auto it = forInsert.begin(); it != forInsert.end(); it++) {	
						triples_for_insert.insert(*it);  //tI_ID_C++; 	
					}
					forInsert.clear();

					// points to the first triple that starts with 's_old' as subject in 'triples_for_insert_2'
					auto ub3 = triples_for_insert_2.upper_bound(DifferentialIndex::SetTriple(s_old,0,0));
					// points to the first triple that starts with '> s_old' as subject in 'triples_for_insert_2'
					auto ub4 = triples_for_insert_2.upper_bound(DifferentialIndex::SetTriple((s_old+1),0,0));
					for (auto it = ub3; it != ub4; it++) {
						unsigned v2, v3;
						v2 = it->predicate; v3 = it->object;
						forMain.push_back(DifferentialIndex::SetTriple(v3,v2,s_old));	
						forInsert.push_back(DifferentialIndex::SetTriple(s_new,v2,v3));
					}
					// Replacement
					if (ub3 != ub4) triples_for_insert_2.erase(ub3,ub4);
					for (auto it = forInsert.begin(); it != forInsert.end(); it++) {
						triples_for_insert_2.insert(*it);  //tI_ID_C++;
					}
					forInsert.clear();

					// 'forSecondary' utilization 
					for (auto it = forSecondary.begin(); it != forSecondary.end(); it++) {
						auto fd = triples_for_insert_2.find(*it);
						if (fd != triples_for_insert_2.end()) {
							triples_for_insert_2.erase(fd); 
							triples_for_insert_2.insert(DifferentialIndex::SetTriple(it->subject,it->predicate,s_new));
							//tI_ID_C++;
						}	
					}
	  				// 'forMain' utilization 
					for (auto it = forMain.begin(); it != forMain.end(); it++) {
						auto fd = triples_for_insert.find(*it);	
						if (fd != triples_for_insert.end()) {		
							triples_for_insert.erase(fd);	
							triples_for_insert.insert(DifferentialIndex::SetTriple(it->subject,it->predicate,s_new));
							//tI_ID_C++;
						}
					}
					// ..........................................................................
					// Update all the 's_old-triples' with the new s value 's_new' in   	
					// both 'triples_for_insert' and 'triples_for_insert_2' triple sets  
					// (we consider safe to not examine the 'predicate' triple places) < END >	
	
					// Update info for Case_1
					notInDB.erase(s_old); notInDB.insert(s_new);
				}

				// CASE_2: 's_old' EXISTS in Database
				else 
				{
					#if UPDATE_LOG
					cout << "CASE_2: COMPLEX" << endl;
					#endif

					//---------------------------------------------------------------------------	
					// STEP_1: ** Find all the s_old-base_triples 'base_triples' **    
					//---------------------------------------------------------------------------
					vector<DifferentialIndex::SetTriple> base_triples;
					// collect all the triples from Database that "s_old" appears as 'subject'
					scan_SPO.collectTriples_S(differentialIndex.getDatabase().getFacts(Database::Order_Subject_Predicate_Object), s_old, base_triples);
					fL_DB_C++;
					// collect all the triples from Database that "s_old" appears as 'object' 
					scan_OPS.collectTriples_O(differentialIndex.getDatabase().getFacts(Database::Order_Object_Predicate_Subject), s_old, base_triples);
					fL_DB_C++;
					//assert(!base_triples.empty());

					//---------------------------------------------------------------------------
					// STEP_2: ** Find all the s_old-clear_triples 'clear_triples' 
					//            = ('base_triples' - 'triples_for_delete') **    
					//---------------------------------------------------------------------------
					vector<DifferentialIndex::SetTriple> clear_triples;
					for (auto it = base_triples.begin(); it != base_triples.end(); it++) {
						auto fd = triples_for_delete.find(*it);
						if (fd == triples_for_delete.end())  // not found
							clear_triples.push_back(*it);	  		
					}
		
					//---------------------------------------------------------------------------	
					// STEP_3: ** Put all the 'clear_triples' to 'triples_for_delete',
					//         'triples_for_insert' and 'triples_for_insert_2' **  
					//---------------------------------------------------------------------------	
					for (auto it = clear_triples.begin(); it != clear_triples.end(); it++) { 
						triples_for_delete.insert(*it);
						triples_for_insert.insert(*it);
						triples_for_insert_2.insert(DifferentialIndex::SetTriple(it->object,it->predicate,it->subject));
					}

					//---------------------------------------------------------------------------	
					// STEP_4: ** Update all the 's_old-triples' with the new s value 's_new' in   	
					//			  both 'triples_for_insert' and 'triples_for_insert_2' triple sets  
					//         (we consider safe to not examine the 'predicate' triple places) **		
					//---------------------------------------------------------------------------	
					vector<DifferentialIndex::SetTriple> forMain, forSecondary;
					vector<DifferentialIndex::SetTriple> forInsert;

					// "Note": for (ub1,ub2) and (ub3,ub4) it holds that if 's_old' does not exist 
					// in 'triples_for_insert' and 'triples_for_insert_2' as subject within a triple, 
					// then ub1 = ub2 and ub3 = ub4, so the respective 'for' will not be executed  

					// points to the first triple that starts with 's_old' as subject in 'triples_for_insert'
					auto ub1 = triples_for_insert.upper_bound(DifferentialIndex::SetTriple(s_old,0,0));
					// points to the first triple that starts with '> s_old' as subject in 'triples_for_insert'
					auto ub2 = triples_for_insert.upper_bound(DifferentialIndex::SetTriple((s_old+1),0,0));
					for (auto it = ub1; it != ub2; it++) {
						unsigned v2, v3;
						v2 = it->predicate; v3 = it->object;
						//#if (!LGD)
						forSecondary.push_back(DifferentialIndex::SetTriple(v3,v2,s_old));
						//#endif	
						forInsert.push_back(DifferentialIndex::SetTriple(s_new,v2,v3));
					}		
					// Replacement
					if (ub1 != ub2) triples_for_insert.erase(ub1,ub2); 
					for (auto it = forInsert.begin(); it != forInsert.end(); it++)	{
						triples_for_insert.insert(*it);  //tI_ID_C++;	
					}
					forInsert.clear();

					// points to the first triple that starts with 's_old' as subject in 'triples_for_insert_2'
					auto ub3 = triples_for_insert_2.upper_bound(DifferentialIndex::SetTriple(s_old,0,0));
					// points to the first triple that starts with '> s_old' as subject in 'triples_for_insert_2'
					auto ub4 = triples_for_insert_2.upper_bound(DifferentialIndex::SetTriple((s_old+1),0,0));
					for (auto it = ub3; it != ub4; it++) {
						unsigned v2, v3;
						v2 = it->predicate; v3 = it->object;
						forMain.push_back(DifferentialIndex::SetTriple(v3,v2,s_old));	
						forInsert.push_back(DifferentialIndex::SetTriple(s_new,v2,v3));
					}
					// Replacement
					if (ub3 != ub4) triples_for_insert_2.erase(ub3,ub4);
					for (auto it = forInsert.begin(); it != forInsert.end(); it++) {
						triples_for_insert_2.insert(*it);  //tI_ID_C++;
					}
					forInsert.clear();

					// 'forSecondary' utilization 
					for (auto it = forSecondary.begin(); it != forSecondary.end(); it++) {
						auto fd = triples_for_insert_2.find(*it);
						if (fd != triples_for_insert_2.end()) {
							triples_for_insert_2.erase(fd); 
							triples_for_insert_2.insert(DifferentialIndex::SetTriple(it->subject,it->predicate,s_new));
							//tI_ID_C++;	
						}	
					}
	  				// 'forMain' utilization 
					for (auto it = forMain.begin(); it != forMain.end(); it++) {
						auto fd = triples_for_insert.find(*it);	
						if (fd != triples_for_insert.end()) {		
							triples_for_insert.erase(fd);	
							triples_for_insert.insert(DifferentialIndex::SetTriple(it->subject,it->predicate,s_new));
							//tI_ID_C++;
						}
					}

					// Update info for Case_2
					notInDB.insert(s_new);
				}   
			}
		}
		// SECOND WAY...
		else if ( (op == true) && (predicate.compare(hasGeo) != 0) ) 
		{
			#if UPDATE_LOG	
			cout << "SECOND WAY..." << endl;
			#endif	
			// We assume that all 'triple PREDICATES' have type URI
			differentialIndex.getDictionary().string2id_updates(predicate, Type::URI, 0, t.predicate);  fL_dict_C++;
			if (t.predicate == ~0u) {
				#if UPDATE_LOG
				cout << "PRED: " << predicate << " NOT in Dict" << endl; 
				#endif
				differentialIndex.getDictionary().updateDictionary(t.predicate, 0, predicate, Type::URI, 0);  fL_dict_C++;
				notInDB.insert(t.predicate);
			} 
			// We assume that all 'triple OBJECTS' have two separate cases
			if (Type::hasSubType(objectType)) {	
				#if UPDATE_LOG	
				cout << "OBJ: " << object << " is NOT a simple string (Type != URI)" << endl;
				#endif
				objectSubType_type = Type::getSubTypeType(objectType);
				differentialIndex.getDictionary().string2id_updates(objectSubType, objectSubType_type, 0, objectSubType_id);  fL_dict_C++;
				if (objectSubType_id == ~0u) {
					#if UPDATE_LOG
					cout << "OBJ_SubType: " << objectSubType << " NOT in Dict" << endl;
					#endif	
					differentialIndex.getDictionary().updateDictionary(objectSubType_id, 0, objectSubType, objectSubType_type, 0);  fL_dict_C++;	
					differentialIndex.getDictionary().updateDictionary(t.object, 0, object, objectType, objectSubType_id);  fL_dict_C++;			
					triples_for_insert.insert(t);
					triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));	
					notInDB.insert(t.object);			
				}	
				else {
					differentialIndex.getDictionary().string2id_updates(object, objectType, objectSubType_id, t.object);  fL_dict_C++;	
					if (t.object == ~0u) {
						#if UPDATE_LOG				
						cout << "OBJ: " << object << " NOT in Dict" << endl;
						#endif
						differentialIndex.getDictionary().updateDictionary(t.object, 0, object, objectType, objectSubType_id);  fL_dict_C++;
						triples_for_insert.insert(t);	
						triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
						notInDB.insert(t.object);
					}
					else {
						#if UPDATE_LOG
						cout << "OBJ: " << object << " EXISTS in Dict" << endl;
						#endif
						if (notInDB.count(t.subject) || notInDB.count(t.predicate) || notInDB.count(t.object)) {
							triples_for_insert.insert(t);		
							triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
						}
						else {
							fL_DB_C++;
							if (!( scan_SPO.findTriple(differentialIndex.getDatabase().getFacts(Database::Order_Subject_Predicate_Object),
							                           t.subject, t.predicate, t.object) )) { 
								triples_for_insert.insert(t);
								triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
							}
						}
					}
				}
			}
			else { 
				#if UPDATE_LOG	
				cout << "OBJ: " << object << " is a simple string (Type = URI or Literal)" << endl;
				#endif 
				differentialIndex.getDictionary().string2id_updates(object, objectType, 0, t.object);  fL_dict_C++;	
				if (t.object == ~0u) {
					#if UPDATE_LOG					
					cout << "OBJ: " << object << " NOT in Dict" << endl; 
					#endif
					differentialIndex.getDictionary().updateDictionary(t.object, 0, object, objectType, 0);  fL_dict_C++;
					triples_for_insert.insert(t);
					triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
					notInDB.insert(t.object);	
				}
				else {
					#if UPDATE_LOG
					cout << "OBJ: " << object << " EXISTS in Dict" << endl;
					#endif
					if (notInDB.count(t.subject) || notInDB.count(t.predicate) || notInDB.count(t.object)) {
						triples_for_insert.insert(t);		
						triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
					}
					else {
						fL_DB_C++;	
						if (!( scan_SPO.findTriple(differentialIndex.getDatabase().getFacts(Database::Order_Subject_Predicate_Object),
						                           t.subject, t.predicate, t.object) )) { 
							triples_for_insert.insert(t);
							triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
						}
					}
				}
			}
		}
		// THIRD WAY... 
		else if ( (op == false) && (predicate.compare(hasGeo) != 0) ) 
		{
			#if UPDATE_LOG
			cout << "THIRD WAY..." << endl;
			#endif
			// We assume that all 'triple PREDICATES' have type URI
			differentialIndex.getDictionary().string2id_updates(predicate, Type::URI, 0, t.predicate);  fL_dict_C++;
			if (t.predicate != ~0u) {
				#if UPDATE_LOG
				cout << "PRED: " << predicate << " EXISTS in Dict" << endl;
				#endif
				// We assume that all 'triple OBJECTS' have two separate cases
				if (Type::hasSubType(objectType)) {	
					#if UPDATE_LOG
					cout << "OBJ: " << object << " is NOT a simple string (Type != URI)" << endl;
					#endif
					objectSubType_type = Type::getSubTypeType(objectType);
					differentialIndex.getDictionary().string2id_updates(objectSubType, objectSubType_type, 0, objectSubType_id);  fL_dict_C++;
					if (objectSubType_id != ~0u) {
						#if UPDATE_LOG
						cout << "OBJ_SubType: " << objectSubType << " EXISTS in Dict" << endl;	
						#endif	
						differentialIndex.getDictionary().string2id_updates(object, objectType, objectSubType_id, t.object);  fL_dict_C++;	
						if (t.object != ~0u) {
							#if UPDATE_LOG				
							cout << "OBJ: " << object << " EXISTS in Dict" << endl;
							#endif
							if (notInDB.count(t.subject) || notInDB.count(t.predicate) || notInDB.count(t.object)) {
								auto fd = triples_for_insert.find(t);
								if (fd != triples_for_insert.end()) { 
									triples_for_insert.erase(fd);  //tI_DEL_C++;
								}
								auto fd_2 = triples_for_insert_2.find(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
								if (fd_2 != triples_for_insert_2.end()) {
									triples_for_insert_2.erase(fd_2);  //tI_DEL_C++;
								}		
							}	
							else {
								fL_DB_C++; 
								if ( scan_SPO.findTriple(differentialIndex.getDatabase().getFacts(Database::Order_Subject_Predicate_Object), 
					                                  t.subject, t.predicate, t.object) ) triples_for_delete.insert(t);
							}
						}					
					}	
				}
				else { 
					#if UPDATE_LOG
					cout << "OBJ: " << object << " is a simple string (Type = URI or Literal)" << endl; 
					#endif
					differentialIndex.getDictionary().string2id_updates(object, objectType, 0, t.object);  fL_dict_C++;	
					if (t.object != ~0u) {
						#if UPDATE_LOG				
						cout << "OBJ: " << object << " EXISTS in Dict" << endl;
						#endif	
						if (notInDB.count(t.subject) || notInDB.count(t.predicate) || notInDB.count(t.object)) {
							auto fd = triples_for_insert.find(t);
							if (fd != triples_for_insert.end()) { 
								triples_for_insert.erase(fd);  //tI_DEL_C++;
 							}
							auto fd_2 = triples_for_insert_2.find(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
							if (fd_2 != triples_for_insert_2.end()) { 
								triples_for_insert_2.erase(fd_2);  //tI_DEL_C++;
							}
						}	
						else {
							fL_DB_C++; 
							if ( scan_SPO.findTriple(differentialIndex.getDatabase().getFacts(Database::Order_Subject_Predicate_Object),
							                         t.subject, t.predicate, t.object) ) triples_for_delete.insert(t);
						}
					}
				}
			} 			
		} 
	}
	// CASE: 's' EXISTS in Dictionary AND 's' is NON-SPATIAL < END >

	// ............................................................................... //

	// CASE: 's' EXISTS in Dictionary AND 's' is SPATIAL < START >
	else if ( (t.subject != ~0u) && ((t.subject % 2) != 0) ) 
	{
		#if UPDATE_LOG
		cout << "** UPDATES (3/3) CASE **" << endl;	
		#endif
		// FIRST WAY...
		if ( (op == true) && (predicate.compare(hasGeo) != 0) ) 
		{
			#if UPDATE_LOG
			cout << "FIRST WAY..." << endl;
			#endif
			// We assume that all 'triple PREDICATES' have type URI
			differentialIndex.getDictionary().string2id_updates(predicate, Type::URI, 0, t.predicate);  fL_dict_C++;
			if (t.predicate == ~0u) {
				#if UPDATE_LOG
				cout << "PRED: " << predicate << " NOT in Dict" << endl; 
				#endif
				differentialIndex.getDictionary().updateDictionary(t.predicate, 0, predicate, Type::URI, 0);  fL_dict_C++;
				notInDB.insert(t.predicate);
			} 
			// We assume that all 'triple OBJECTS' have two separate cases
			if (Type::hasSubType(objectType)) {
				#if UPDATE_LOG	
				cout << "OBJ: " << object << " is NOT a simple string (Type != URI)" << endl;
				#endif
				objectSubType_type = Type::getSubTypeType(objectType);
				differentialIndex.getDictionary().string2id_updates(objectSubType, objectSubType_type, 0, objectSubType_id);  fL_dict_C++;
				if (objectSubType_id == ~0u) {
					#if UPDATE_LOG
					cout << "OBJ_SubType: " << objectSubType << " NOT in Dict" << endl;	
					#endif
					differentialIndex.getDictionary().updateDictionary(objectSubType_id, 0, objectSubType, objectSubType_type, 0);  fL_dict_C++;	
					differentialIndex.getDictionary().updateDictionary(t.object, 0, object, objectType, objectSubType_id);  fL_dict_C++;			
					triples_for_insert.insert(t);	
					triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
					notInDB.insert(t.object);			
				}	
				else {
					differentialIndex.getDictionary().string2id_updates(object, objectType, objectSubType_id, t.object);  fL_dict_C++;	
					if (t.object == ~0u) {
						#if UPDATE_LOG				
						cout << "OBJ: " << object << " NOT in Dict" << endl;
						#endif
						differentialIndex.getDictionary().updateDictionary(t.object, 0, object, objectType, objectSubType_id);  fL_dict_C++;
						triples_for_insert.insert(t);	
						triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
						notInDB.insert(t.object);
					}
					else {
						#if UPDATE_LOG
						cout << "OBJ: " << object << " EXISTS in Dict" << endl;
						#endif
						if (notInDB.count(t.subject) || notInDB.count(t.predicate) || notInDB.count(t.object)) {
							triples_for_insert.insert(t);		
							triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
						} 
						else {
							fL_DB_C++; 
							if (!( scan_SPO.findTriple(differentialIndex.getDatabase().getFacts(Database::Order_Subject_Predicate_Object),
							                           t.subject, t.predicate, t.object) )) {
								triples_for_insert.insert(t);
								triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
							}
						}			
					}
				}
			}
			else { 
				#if UPDATE_LOG
				cout << "OBJ: " << object << " is a simple string (Type = URI or Literal)" << endl;
				#endif 
				differentialIndex.getDictionary().string2id_updates(object, objectType, 0, t.object);  fL_dict_C++;	
				if (t.object == ~0u) {
					#if UPDATE_LOG				
					cout << "OBJ: " << object << " NOT in Dict" << endl; 
					#endif
					differentialIndex.getDictionary().updateDictionary(t.object, 0, object, objectType, 0);  fL_dict_C++; 
					triples_for_insert.insert(t);
					triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
					notInDB.insert(t.object);
				}
				else {
					#if UPDATE_LOG
					cout << "OBJ: " << object << " EXISTS in Dict" << endl;
					#endif
					if (notInDB.count(t.subject) || notInDB.count(t.predicate) || notInDB.count(t.object)) { 
						triples_for_insert.insert(t);
						triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
					}
					else {
						fL_DB_C++; 
						if (!( scan_SPO.findTriple(differentialIndex.getDatabase().getFacts(Database::Order_Subject_Predicate_Object),
						                           t.subject, t.predicate, t.object) )) { 
							triples_for_insert.insert(t);
							triples_for_insert_2.insert(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
						}
					}
				}
			}
		}
		// SECOND WAY...
		else if ( (op == false) && (predicate.compare(hasGeo) == 0) )			
		{
			#if UPDATE_LOG
			cout << "SECOND WAY..." << endl;
			#endif
			// We assume that all 'triple PREDICATES' have type URI
			differentialIndex.getDictionary().string2id_updates(predicate, Type::URI, 0, t.predicate);  fL_dict_C++;
			if (t.predicate != ~0u) {
				#if UPDATE_LOG
				cout << "PRED: " << predicate << " EXISTS in Dict" << endl;
				#endif
				objectSubType_type = Type::getSubTypeType(objectType);
				differentialIndex.getDictionary().string2id_updates(objectSubType, objectSubType_type, 0, objectSubType_id);  fL_dict_C++;
				if (objectSubType_id != ~0u) {
					#if UPDATE_LOG
					cout << "OBJ_SubType: " << objectSubType << " EXISTS in Dict" << endl;	
					#endif
					differentialIndex.getDictionary().string2id_updates(object, objectType, objectSubType_id, t.object);  fL_dict_C++;	
					if (t.object != ~0u) 
					{
						#if UPDATE_LOG					
						cout << "OBJ: " << object << " EXISTS in Dict" << endl;
						#endif
						// Flag that checks if we must skip or not the resolved triple 't'
						bool skip = 1;

						if (notInDB.count(t.subject)) {
							auto fd = triples_for_insert.find(t);
							if (fd != triples_for_insert.end()) {
								triples_for_insert.erase(fd);  skip = 0;  //tI_DEL_C++; 	
							}  								
						}
						else {
							fL_DB_C++; 
							if ( scan_SPO.findTriple(differentialIndex.getDatabase().getFacts(Database::Order_Subject_Predicate_Object),
							                         t.subject, t.predicate, t.object) ) { 
								triples_for_delete.insert(t);  skip = 0;
							}
						}
			
						if (!skip)
						{
						// each entry is a changed ID pair (old, new)
						vector<pair<unsigned,unsigned>> changes;  

						unsigned s_old = t.subject; 
						unsigned s_new = differentialIndex.getDictionary().maxEvenID + 2; 
						differentialIndex.getDictionary().maxEvenID = s_new; 	
	
						changes.push_back(make_pair(s_old, s_new));
						differentialIndex.getDictionary().updateDictForIDs(subject, Type::URI, 0, s_old, s_new);  fL_dict_C++;

						// get the upper level of 's_old' level
						unsigned up_lvl;  
						unsigned level = ((s_old & LVMASK) >> 1);  // 's_old' level
						if (level == 0) up_lvl = 0;
						else up_lvl = level - 1;  

						differentialIndex.getDictionary().reEncodeGridAfterDelete(
							s_old, up_lvl, 1, differentialIndex.getDatabase().grid, 
							differentialIndex.getDatabase().getFacts(Database::Order_Subject_Predicate_Object),						
							t.predicate, notInDB, changes, fL_DB_C, fL_dict_C
						);

						// ....
						// Right now 'changes' must have been filled OK.
						// Apply the change ID code in updates(2/3) first way
						// ....

						for (auto it_C = changes.begin(); it_C != changes.end(); it_C++)
						{
							unsigned s_old = it_C->first;  unsigned s_new = it_C->second;
							ReEncode_C++;

							// CASE_1: "s_old" NOT in Database
							if (notInDB.count(s_old))
							{
								#if UPDATE_LOG
								cout << "CASE_1: SIMPLE" << endl;
								#endif	

								// Update all the 's_old-triples' with the new s value 's_new' in   	
								// both 'triples_for_insert' and 'triples_for_insert_2' triple sets  
								// (we consider safe to not examine the 'predicate' triple places) < START >	 
								// ..........................................................................
								vector<DifferentialIndex::SetTriple> forMain, forSecondary;
								vector<DifferentialIndex::SetTriple> forInsert;

								// "Note": for (ub1,ub2) and (ub3,ub4) it holds that if 's_old' does not exist 
								// in 'triples_for_insert' and 'triples_for_insert_2' as subject within a triple, 
								// then ub1 = ub2 and ub3 = ub4, so the respective 'for' will not be executed  

								// points to the first triple that starts with 's_old' as subject in 'triples_for_insert'
								auto ub1 = triples_for_insert.upper_bound(DifferentialIndex::SetTriple(s_old,0,0));
								// points to the first triple that starts with '> s_old' as subject in 'triples_for_insert'
								auto ub2 = triples_for_insert.upper_bound(DifferentialIndex::SetTriple((s_old+1),0,0));
								for (auto it = ub1; it != ub2; it++) {
									unsigned v2, v3;
									v2 = it->predicate; v3 = it->object;
									//#if (!LGD)
									forSecondary.push_back(DifferentialIndex::SetTriple(v3,v2,s_old));
									//#endif	
									forInsert.push_back(DifferentialIndex::SetTriple(s_new,v2,v3));
								}		
								// Replacement
								if (ub1 != ub2) triples_for_insert.erase(ub1,ub2); 
								for (auto it = forInsert.begin(); it != forInsert.end(); it++) {	
									triples_for_insert.insert(*it);  //tI_ID_C++;	
								}
								forInsert.clear();

								// points to the first triple that starts with 's_old' as subject in 'triples_for_insert_2'
								auto ub3 = triples_for_insert_2.upper_bound(DifferentialIndex::SetTriple(s_old,0,0));
								// points to the first triple that starts with '> s_old' as subject in 'triples_for_insert_2'
								auto ub4 = triples_for_insert_2.upper_bound(DifferentialIndex::SetTriple((s_old+1),0,0));
								for (auto it = ub3; it != ub4; it++) {
									unsigned v2, v3;
									v2 = it->predicate; v3 = it->object;
									forMain.push_back(DifferentialIndex::SetTriple(v3,v2,s_old));	
									forInsert.push_back(DifferentialIndex::SetTriple(s_new,v2,v3));
								}
								// Replacement
								if (ub3 != ub4) triples_for_insert_2.erase(ub3,ub4);
								for (auto it = forInsert.begin(); it != forInsert.end(); it++) {
									triples_for_insert_2.insert(*it);  //tI_ID_C++;
								}
								forInsert.clear();

								// 'forSecondary' utilization 
								for (auto it = forSecondary.begin(); it != forSecondary.end(); it++) {
									auto fd = triples_for_insert_2.find(*it);
									if (fd != triples_for_insert_2.end()) {
										triples_for_insert_2.erase(fd); 
										triples_for_insert_2.insert(DifferentialIndex::SetTriple(it->subject,it->predicate,s_new));
										//tI_ID_C++;
									}	
								}
	  							// 'forMain' utilization 
								for (auto it = forMain.begin(); it != forMain.end(); it++) {
									auto fd = triples_for_insert.find(*it);	
									if (fd != triples_for_insert.end()) {		
										triples_for_insert.erase(fd);	
										triples_for_insert.insert(DifferentialIndex::SetTriple(it->subject,it->predicate,s_new));
										//tI_ID_C++;			
									}
								}
								// ..........................................................................
								// Update all the 's_old-triples' with the new s value 's_new' in   	
								// both 'triples_for_insert' and 'triples_for_insert_2' triple sets  
								// (we consider safe to not examine the 'predicate' triple places) < END >	
	
								// Update info for Case_1
								notInDB.erase(s_old); notInDB.insert(s_new);
							}
	
							// CASE_2: 's_old' EXISTS in Database
							else 
							{
								#if UPDATE_LOG
								cout << "CASE_2: COMPLEX" << endl;
								#endif

								//---------------------------------------------------------------------------	
								// STEP_1: ** Find all the s_old-base_triples 'base_triples' **    
								//---------------------------------------------------------------------------
								vector<DifferentialIndex::SetTriple> base_triples;
								// collect all the triples from Database that "s_old" appears as 'subject'
								scan_SPO.collectTriples_S(differentialIndex.getDatabase().getFacts(Database::Order_Subject_Predicate_Object), s_old, base_triples);
								fL_DB_C++;	
								// collect all the triples from Database that "s_old" appears as 'object' 
								scan_OPS.collectTriples_O(differentialIndex.getDatabase().getFacts(Database::Order_Object_Predicate_Subject), s_old, base_triples);
								fL_DB_C++;
								//assert(!base_triples.empty());

								//---------------------------------------------------------------------------
								// STEP_2: ** Find all the s_old-clear_triples 'clear_triples' 
								//            = ('base_triples' - 'triples_for_delete') **    
								//---------------------------------------------------------------------------
								vector<DifferentialIndex::SetTriple> clear_triples;
								for (auto it = base_triples.begin(); it != base_triples.end(); it++) {
									auto fd = triples_for_delete.find(*it);
									if (fd == triples_for_delete.end())  // not found
										clear_triples.push_back(*it);	  		
								}
		
								//---------------------------------------------------------------------------	
								// STEP_3: ** Put all the 'clear_triples' to 'triples_for_delete',
								//         'triples_for_insert' and 'triples_for_insert_2' **  
								//---------------------------------------------------------------------------	
								for (auto it = clear_triples.begin(); it != clear_triples.end(); it++) { 
									triples_for_delete.insert(*it);
									triples_for_insert.insert(*it);
									triples_for_insert_2.insert(DifferentialIndex::SetTriple(it->object,it->predicate,it->subject));
								}

								//---------------------------------------------------------------------------	
								// STEP_4: ** Update all the 's_old-triples' with the new s value 's_new' in   	
								//			  both 'triples_for_insert' and 'triples_for_insert_2' triple sets  
								//         (we consider safe to not examine the 'predicate' triple places) **		
								//---------------------------------------------------------------------------	
								vector<DifferentialIndex::SetTriple> forMain, forSecondary;
								vector<DifferentialIndex::SetTriple> forInsert;

								// "Note": for (ub1,ub2) and (ub3,ub4) it holds that if 's_old' does not exist 
								// in 'triples_for_insert' and 'triples_for_insert_2' as subject within a triple, 
								// then ub1 = ub2 and ub3 = ub4, so the respective 'for' will not be executed  

								// points to the first triple that starts with 's_old' as subject in 'triples_for_insert'
								auto ub1 = triples_for_insert.upper_bound(DifferentialIndex::SetTriple(s_old,0,0));
								// points to the first triple that starts with '> s_old' as subject in 'triples_for_insert'
								auto ub2 = triples_for_insert.upper_bound(DifferentialIndex::SetTriple((s_old+1),0,0));
								for (auto it = ub1; it != ub2; it++) {
									unsigned v2, v3;
									v2 = it->predicate; v3 = it->object;
									//#if (!LGD)
									forSecondary.push_back(DifferentialIndex::SetTriple(v3,v2,s_old));
									//#endif	
									forInsert.push_back(DifferentialIndex::SetTriple(s_new,v2,v3));
								}		
								// Replacement
								if (ub1 != ub2) triples_for_insert.erase(ub1,ub2); 
								for (auto it = forInsert.begin(); it != forInsert.end(); it++) {	
									triples_for_insert.insert(*it);  //tI_ID_C++;	
								}
								forInsert.clear();

								// points to the first triple that starts with 's_old' as subject in 'triples_for_insert_2'
								auto ub3 = triples_for_insert_2.upper_bound(DifferentialIndex::SetTriple(s_old,0,0));
								// points to the first triple that starts with '> s_old' as subject in 'triples_for_insert_2'
								auto ub4 = triples_for_insert_2.upper_bound(DifferentialIndex::SetTriple((s_old+1),0,0));
								for (auto it = ub3; it != ub4; it++) {
									unsigned v2, v3;
									v2 = it->predicate; v3 = it->object;
									forMain.push_back(DifferentialIndex::SetTriple(v3,v2,s_old));	
									forInsert.push_back(DifferentialIndex::SetTriple(s_new,v2,v3));
								}
								// Replacement
								if (ub3 != ub4) triples_for_insert_2.erase(ub3,ub4);
								for (auto it = forInsert.begin(); it != forInsert.end(); it++) {
									triples_for_insert_2.insert(*it);  //tI_ID_C++;
								}
								forInsert.clear();

								// 'forSecondary' utilization 
								for (auto it = forSecondary.begin(); it != forSecondary.end(); it++) {
									auto fd = triples_for_insert_2.find(*it);
									if (fd != triples_for_insert_2.end()) {
										triples_for_insert_2.erase(fd); 
										triples_for_insert_2.insert(DifferentialIndex::SetTriple(it->subject,it->predicate,s_new));
										//tI_ID_C++;
									}	
								}
	  							// 'forMain' utilization 
								for (auto it = forMain.begin(); it != forMain.end(); it++) {
									auto fd = triples_for_insert.find(*it);	
									if (fd != triples_for_insert.end()) {		
										triples_for_insert.erase(fd);	
										triples_for_insert.insert(DifferentialIndex::SetTriple(it->subject,it->predicate,s_new));
										//tI_ID_C++;	
									}
								}

								// Update info for Case_2
								notInDB.insert(s_new);
							}		  
						}
						}
					}
				}
			}
		}
		// THIRD WAY...
		else if ( (op == false) && (predicate.compare(hasGeo) != 0) )
		{
			#if UPDATE_LOG
			cout << "THIRD WAY..." << endl;
			#endif	
			// We assume that all 'triple PREDICATES' have type URI
			differentialIndex.getDictionary().string2id_updates(predicate, Type::URI, 0, t.predicate);  fL_dict_C++;
			if (t.predicate != ~0u) {
				#if UPDATE_LOG
				cout << "PRED: " << predicate << " EXISTS in Dict" << endl;
				#endif
				// We assume that all 'triple OBJECTS' have two separate cases
				if (Type::hasSubType(objectType)) {
					#if UPDATE_LOG		
					cout << "OBJ: " << object << " is NOT a simple string (Type != URI)" << endl;
					#endif
					objectSubType_type = Type::getSubTypeType(objectType);
					differentialIndex.getDictionary().string2id_updates(objectSubType, objectSubType_type, 0, objectSubType_id);  fL_dict_C++;
					if (objectSubType_id != ~0u) {
						#if UPDATE_LOG
						cout << "OBJ_SubType: " << objectSubType << " EXISTS in Dict" << endl;	
						#endif
						differentialIndex.getDictionary().string2id_updates(object, objectType, objectSubType_id, t.object);  fL_dict_C++;	
						if (t.object != ~0u) {
							#if UPDATE_LOG					
							cout << "OBJ: " << object << " EXISTS in Dict" << endl;
							#endif
							if (notInDB.count(t.subject) || notInDB.count(t.predicate) || notInDB.count(t.object)) {	
								auto fd = triples_for_insert.find(t);
								if (fd != triples_for_insert.end()) { 
									triples_for_insert.erase(fd);  //tI_DEL_C++;
								}
								auto fd_2 = triples_for_insert_2.find(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
								if (fd_2 != triples_for_insert_2.end()) { 
									triples_for_insert_2.erase(fd_2);  //tI_DEL_C++;
								}
							}	
							else {
								fL_DB_C++;
								if ( scan_SPO.findTriple(differentialIndex.getDatabase().getFacts(Database::Order_Subject_Predicate_Object), 
					                                  t.subject, t.predicate, t.object) ) triples_for_delete.insert(t);
							}
						}					
					}	
				}
				else { 
					#if UPDATE_LOG
					cout << "OBJ: " << object << " is a simple string (Type = URI or Literal)" << endl; 
					#endif
					differentialIndex.getDictionary().string2id_updates(object, objectType, 0, t.object);  fL_dict_C++;	
					if (t.object != ~0u) {				
						#if UPDATE_LOG	
						cout << "OBJ: " << object << " EXISTS in Dict" << endl;
						#endif
						if (notInDB.count(t.subject) || notInDB.count(t.predicate) || notInDB.count(t.object)) {
							auto fd = triples_for_insert.find(t);
							if (fd != triples_for_insert.end()) { 
								triples_for_insert.erase(fd);  //tI_DEL_C++; 
							}
							auto fd_2 = triples_for_insert_2.find(DifferentialIndex::SetTriple(t.object,t.predicate,t.subject));
							if (fd_2 != triples_for_insert_2.end()) { 
								triples_for_insert_2.erase(fd_2);  //tI_DEL_C++;
							} 	
						}	
						else { 
							fL_DB_C++;
							if ( scan_SPO.findTriple(differentialIndex.getDatabase().getFacts(Database::Order_Subject_Predicate_Object),
							                         t.subject, t.predicate, t.object) ) triples_for_delete.insert(t);
						}	
					}
				}
			} 
		}
	}
	// CASE: 's' EXISTS in Dictionary AND 's' is SPATIAL < END >
#endif
}
//---------------------------------------------------------------------------
namespace {
//---------------------------------------------------------------------------
/// Sort by subject
struct SortBySubject { bool operator()(const DifferentialIndex::Triple& a,const DifferentialIndex::Triple& b) const { return a.subject<b.subject; } };
/// Sort by predicate
struct SortByPredicate { bool operator()(const DifferentialIndex::Triple& a,const DifferentialIndex::Triple& b) const { return a.predicate<b.predicate; } };
/// Sort by object
struct SortByObject { bool operator()(const DifferentialIndex::Triple& a,const DifferentialIndex::Triple& b) const { return a.object<b.object; } };
//---------------------------------------------------------------------------
static double boxVolume(const PredicateLockManager::Box& b)
   // Compute the volume of a box
{
   return (static_cast<double>(b.subjectMax-b.subjectMin)+1.0)*
          (static_cast<double>(b.predicateMax-b.predicateMin)+1.0)*
          (static_cast<double>(b.objectMax-b.objectMin)+1.0);
}
//---------------------------------------------------------------------------
static double computeSplitVolume(vector<DifferentialIndex::Triple>::iterator iter,vector<DifferentialIndex::Triple>::iterator limit,unsigned split,unsigned slot)
   // Compute the volumes after a split
{
   PredicateLockManager::Box leftBox(0,0,0,0,0,0),rightBox(0,0,0,0,0,0);
   bool firstLeft=true,firstRight=true;

   for (;iter!=limit;++iter) {
      // Determine the partition
      bool left=false;
      switch (slot) {
         case 0: left=(*iter).subject<=split; break;
         case 1: left=(*iter).predicate<=split; break;
         case 2: left=(*iter).object<=split; break;
      }
      // Remember
      if (left) {
         if (((*iter).subject<leftBox.subjectMin)||(firstLeft))
            leftBox.subjectMin=(*iter).subject;
         if (((*iter).subject>leftBox.subjectMax)||(firstLeft))
            leftBox.subjectMax=(*iter).subject;
         if (((*iter).predicate<leftBox.predicateMin)||(firstLeft))
            leftBox.predicateMin=(*iter).predicate;
         if (((*iter).predicate>leftBox.predicateMax)||(firstLeft))
            leftBox.predicateMax=(*iter).predicate;
         if (((*iter).object<leftBox.objectMin)||(firstLeft))
            leftBox.objectMin=(*iter).object;
         if (((*iter).object>leftBox.objectMax)||(firstLeft))
            leftBox.objectMax=(*iter).object;
         firstLeft=false;
      } else {
         if (((*iter).subject<rightBox.subjectMin)||(firstRight))
            rightBox.subjectMin=(*iter).subject;
         if (((*iter).subject>rightBox.subjectMax)||(firstRight))
            rightBox.subjectMax=(*iter).subject;
         if (((*iter).predicate<rightBox.predicateMin)||(firstRight))
            rightBox.predicateMin=(*iter).predicate;
         if (((*iter).predicate>rightBox.predicateMax)||(firstRight))
            rightBox.predicateMax=(*iter).predicate;
         if (((*iter).object<rightBox.objectMin)||(firstRight))
            rightBox.objectMin=(*iter).object;
         if (((*iter).object>rightBox.objectMax)||(firstRight))
            rightBox.objectMax=(*iter).object;
         firstRight=false;
      }
   }
   return boxVolume(leftBox)+boxVolume(rightBox);
}
//---------------------------------------------------------------------------
static vector<DifferentialIndex::Triple>::iterator computeSplit(vector<DifferentialIndex::Triple>::iterator iter,vector<DifferentialIndex::Triple>::iterator limit,unsigned splitValue,unsigned slot,PredicateLockManager::Box& leftBox,PredicateLockManager::Box& rightBox)
   // Compute the volumes after a split
{
   bool firstLeft=true,firstRight=true;

   vector<DifferentialIndex::Triple>::iterator split=iter;
   for (;iter!=limit;++iter) {
      // Determine the partition
      bool left=false;
      switch (slot) {
         case 0: left=(*iter).subject<=splitValue; break;
         case 1: left=(*iter).predicate<=splitValue; break;
         case 2: left=(*iter).object<=splitValue; break;
      }
      // Remember
      if (left) {
         if (((*iter).subject<leftBox.subjectMin)||(firstLeft))
            leftBox.subjectMin=(*iter).subject;
         if (((*iter).subject>leftBox.subjectMax)||(firstLeft))
            leftBox.subjectMax=(*iter).subject;
         if (((*iter).predicate<leftBox.predicateMin)||(firstLeft))
            leftBox.predicateMin=(*iter).predicate;
         if (((*iter).predicate>leftBox.predicateMax)||(firstLeft))
            leftBox.predicateMax=(*iter).predicate;
         if (((*iter).object<leftBox.objectMin)||(firstLeft))
            leftBox.objectMin=(*iter).object;
         if (((*iter).object>leftBox.objectMax)||(firstLeft))
            leftBox.objectMax=(*iter).object;
         firstLeft=false;
         swap(*iter,*split);
         ++split;
      } else {
         if (((*iter).subject<rightBox.subjectMin)||(firstRight))
            rightBox.subjectMin=(*iter).subject;
         if (((*iter).subject>rightBox.subjectMax)||(firstRight))
            rightBox.subjectMax=(*iter).subject;
         if (((*iter).predicate<rightBox.predicateMin)||(firstRight))
            rightBox.predicateMin=(*iter).predicate;
         if (((*iter).predicate>rightBox.predicateMax)||(firstRight))
            rightBox.predicateMax=(*iter).predicate;
         if (((*iter).object<rightBox.objectMin)||(firstRight))
            rightBox.objectMin=(*iter).object;
         if (((*iter).object>rightBox.objectMax)||(firstRight))
            rightBox.objectMax=(*iter).object;
         firstRight=false;
      }
   }
   return split;
}
//---------------------------------------------------------------------------
}
//---------------------------------------------------------------------------
void BulkOperation::buildCover(unsigned maxSize,vector<PredicateLockManager::Box>& boxes)
{
   boxes.clear();

   // Trivial?
   if (triples.size()<=maxSize) {
      for (vector<DifferentialIndex::Triple>::iterator iter=triples.begin(),limit=triples.end();iter!=limit;++iter)
         boxes.push_back(PredicateLockManager::Box((*iter).subject,(*iter).subject,(*iter).predicate,(*iter).predicate,(*iter).object,(*iter).object));
      return;
   }

   // Build one big initial box
   vector<vector<DifferentialIndex::Triple>::iterator> lowerBounds,upperBounds;
   lowerBounds.push_back(triples.begin());
   upperBounds.push_back(triples.end());
   boxes.push_back(PredicateLockManager::Box(0,~0u,0,~0u,0,~0u));
   while (boxes.size()<maxSize) {
      // Find the largest box
      unsigned largestBox=0; double largestVolume=boxVolume(boxes[0]);
      for (unsigned index=1,limit=boxes.size();index<limit;index++) {
         double v=boxVolume(boxes[index]);
         if (v>largestVolume) {
            largestBox=index;
            largestVolume=v;
         }
      }
      // The the best split candidate
      unsigned lastValue,bestGap;

      sort(lowerBounds[largestBox],upperBounds[largestBox],SortBySubject());
      unsigned bestSplitSubject=(*lowerBounds[largestBox]).subject;
      lastValue=bestSplitSubject; bestGap=0;
      for (vector<DifferentialIndex::Triple>::iterator iter=lowerBounds[largestBox],limit=upperBounds[largestBox];iter!=limit;++iter) {
         unsigned gap=(*iter).subject-lastValue;
         lastValue=(*iter).subject;
         if (gap>bestGap) {
            bestSplitSubject=lastValue;
            bestGap=gap;
         }
      }
      double subjectVolume=computeSplitVolume(lowerBounds[largestBox],upperBounds[largestBox],bestSplitSubject,0);

      sort(lowerBounds[largestBox],upperBounds[largestBox],SortByPredicate());
      unsigned bestSplitPredicate=(*lowerBounds[largestBox]).predicate;
      lastValue=bestSplitPredicate; bestGap=0;
      for (vector<DifferentialIndex::Triple>::iterator iter=lowerBounds[largestBox],limit=upperBounds[largestBox];iter!=limit;++iter) {
         unsigned gap=(*iter).predicate-lastValue;
         lastValue=(*iter).predicate;
         if (gap>bestGap) {
            bestSplitPredicate=lastValue;
            bestGap=gap;
         }
      }
      double predicateVolume=computeSplitVolume(lowerBounds[largestBox],upperBounds[largestBox],bestSplitPredicate,1);

      sort(lowerBounds[largestBox],upperBounds[largestBox],SortByObject());
      unsigned bestSplitObject=(*lowerBounds[largestBox]).object;
      lastValue=bestSplitObject; bestGap=0;
      for (vector<DifferentialIndex::Triple>::iterator iter=lowerBounds[largestBox],limit=upperBounds[largestBox];iter!=limit;++iter) {
         unsigned gap=(*iter).object-lastValue;
         lastValue=(*iter).object;
         if (gap>bestGap) {
            bestSplitObject=lastValue;
            bestGap=gap;
         }
      }
      double objectVolume=computeSplitVolume(lowerBounds[largestBox],upperBounds[largestBox],bestSplitSubject,2);

      // Perform the split
      PredicateLockManager::Box leftBox(0,0,0,0,0,0),rightBox(0,0,0,0,0,0);
      vector<DifferentialIndex::Triple>::iterator split;
      if ((predicateVolume<=subjectVolume)&&(predicateVolume<=objectVolume)) {
         split=computeSplit(lowerBounds[largestBox],upperBounds[largestBox],bestSplitPredicate,1,leftBox,rightBox);
      } else if ((subjectVolume<=predicateVolume)&&(subjectVolume<=objectVolume)) {
         split=computeSplit(lowerBounds[largestBox],upperBounds[largestBox],bestSplitSubject,0,leftBox,rightBox);
      } else {
         split=computeSplit(lowerBounds[largestBox],upperBounds[largestBox],bestSplitObject,2,leftBox,rightBox);
      }

      // Update the bounds
      lowerBounds.push_back(split);
      upperBounds.push_back(upperBounds[largestBox]);
      upperBounds[largestBox]=split;
      boxes[largestBox]=leftBox;
      boxes.push_back(rightBox);
   }
}
//---------------------------------------------------------------------------
void BulkOperation::commit()
   // Commit
{
   // Resolve all temporary ids
   vector<unsigned> realIds;
   differentialIndex.mapLiterals(id2string, realIds);
   unsigned tempStart = (~0u) - realIds.size();
   for (vector<DifferentialIndex::Triple>::iterator iter = triples.begin(), limit = triples.end(); iter != limit; ++iter) {
      if ((*iter).subject > tempStart) 
			(*iter).subject = realIds[(~0u) - (*iter).subject];
      if ((*iter).predicate > tempStart) 
			(*iter).predicate = realIds[(~0u) - (*iter).predicate];
      if ((*iter).object > tempStart) 
			(*iter).object = realIds[(~0u) - (*iter).object];
   }
   realIds.clear();

   // Load the triples
   differentialIndex.load(triples, deleteMarker);

   // And release
   id2string.clear();
   string2id.clear();
   triples.clear();
}
//---------------------------------------------------------------------------
void BulkOperation::commit_updates(fstream &f2)
   // Commit (updates version)
{
	// Load the 'DI' with the triples for 'delete' from the database indexes
	if (!triples_for_delete.empty()) {
		deleteMarker = true;  tD_C = 6 * triples_for_delete.size();
		differentialIndex.load_updates(triples_for_delete, deleteMarker);
		//#if UPDATE_LOG
		#if 0
		// Print the 'triples_for_delete' 
		cout << endl << "BO::tD...." << endl;
		for (auto it = triples_for_delete.begin(); it != triples_for_delete.end(); it++)
		   cout << "SUB: " << it->subject << ", PRED: " << it->predicate << ", OBJ: " << it->object << endl;
		#endif
		// Release the 'triples_for_delete' 
		triples_for_delete.clear();
	}

	// Load the 'DI' with the triples for 'insert' to the database indexes
	if (!triples_for_insert.empty()) {
		deleteMarker = false;  tI_C = 6 * triples_for_insert.size();
		differentialIndex.load_updates(triples_for_insert, deleteMarker);
		//#if UPDATE_LOG	
		#if 0
		// Print the 'triples_for_insert'
		cout << endl << "BO::tI...." << endl;
		for (auto it = triples_for_insert.begin(); it != triples_for_insert.end(); it++)
		   cout << "SUB: " << it->subject << ", PRED: " << it->predicate << ", OBJ: " << it->object << endl;
		cout << "TOTAL: " << triples_for_insert.size() << endl;
		#endif
		#if 0		
		// Print the 'triples_for_insert_2'
		cout << endl << "BO::tI_2...." << endl;
		for (auto it = triples_for_insert_2.begin(); it != triples_for_insert_2.end(); it++)
		   cout << "SUB: " << it->subject << ", PRED: " << it->predicate << ", OBJ: " << it->object << endl;
		#endif	
		// Release the 'triples_for_insert'
		triples_for_insert.clear();
		// Release the 'triples_for_insert_2'
		triples_for_insert_2.clear();	
	}

	// Clear Other Batch Structures
	notInDB.clear();

	// Keep stats <Batch_ID, fL_DB_C, fL_dict_C, tI_C, tD_C, ReEncode_C>
	f2 << cntB << " " << fL_DB_C << " " << fL_dict_C << " " << tI_C << " " << tD_C << " " << ReEncode_C;
}
//---------------------------------------------------------------------------
void BulkOperation::abort()
   // Abort
{
   id2string.clear();
   string2id.clear();
   triples.clear();
}
//---------------------------------------------------------------------------
