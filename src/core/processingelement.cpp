/*************************************************************************************/
/*      Copyright 2009 Barcelona Supercomputing Center                               */
/*                                                                                   */
/*      This file is part of the NANOS++ library.                                    */
/*                                                                                   */
/*      NANOS++ is free software: you can redistribute it and/or modify              */
/*      it under the terms of the GNU Lesser General Public License as published by  */
/*      the Free Software Foundation, either version 3 of the License, or            */
/*      (at your option) any later version.                                          */
/*                                                                                   */
/*      NANOS++ is distributed in the hope that it will be useful,                   */
/*      but WITHOUT ANY WARRANTY; without even the implied warranty of               */
/*      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                */
/*      GNU Lesser General Public License for more details.                          */
/*                                                                                   */
/*      You should have received a copy of the GNU Lesser General Public License     */
/*      along with NANOS++.  If not, see <http://www.gnu.org/licenses/>.             */
/*************************************************************************************/

#include <string.h>
#include "processingelement.hpp"
#include "basethread.hpp"
#include "debug.hpp"
#include "schedule.hpp"
#include "copydata.hpp"
#include "system.hpp"
#include "instrumentation.hpp"
#include "directory.hpp"

using namespace nanos;

void ProcessingElement::copyDataIn( WorkDescriptor &work )
{
   Directory *dir = work.getParent()->getDirectory(true);
   if ( dir != NULL ) {
      CopyData *copies = work.getCopies();
      for ( unsigned int i = 0; i < work.getNumCopies(); i++ ) {
         CopyData & cd = copies[i];
         if ( !cd.isPrivate() ) {
              dir->registerAccess( cd.getAddress(), cd.getSize(), cd.isInput(), cd.isOutput() );
         }
      }
   }
}

void ProcessingElement::copyDataOut( WorkDescriptor &work )
{
   Directory *dir = work.getParent()->getDirectory(false);
   if ( dir != NULL ) {
      CopyData *copies = work.getCopies();
      for ( unsigned int i = 0; i < work.getNumCopies(); i++ ) {
         CopyData & cd = copies[i];
         if ( !cd.isPrivate() ) {
            dir->unRegisterAccess( cd.getAddress(), cd.isOutput(), work.getDirectory(false) );
            if ( cd.isOutput() ) {
               Directory *sons = work.getDirectory(false);
               if ( sons!=NULL ) {
                  dir->updateCurrentDirectory( cd.getAddress(), *sons );
               }
            }
         }
      }
   }
}

void ProcessingElement::waitInputs( WorkDescriptor &work )
{
   Directory *dir = work.getParent()->getDirectory(false);
   if ( dir != NULL ) {
      CopyData *copies = work.getCopies();
      for ( unsigned int i = 0; i < work.getNumCopies(); i++ ) {
         CopyData & cd = copies[i];
         if ( !cd.isPrivate() && cd.isInput() ) {
              dir->waitInput( cd.getAddress(), cd.isOutput() );
         }
      }
   }
}

BaseThread& ProcessingElement::startWorker ( )
{
   WD & worker = getWorkerWD();

   NANOS_INSTRUMENT (sys.getInstrumentation()->raiseOpenPtPEvent ( NANOS_WD_DOMAIN, (nanos_event_id_t) worker.getId(), 0, 0 ); )
   NANOS_INSTRUMENT (InstrumentationContextData *icd = worker.getInstrumentationContextData() );
   NANOS_INSTRUMENT (icd->setStartingWD(true) );

   return startThread( worker );
}

BaseThread & ProcessingElement::startThread ( WD &work )
{
   BaseThread &thread = createThread( work );

   thread.start();

   _threads.push_back( &thread );

   return thread;
}

BaseThread & ProcessingElement::associateThisThread ( bool untieMain )
{
   WD & worker = getMasterWD();
   NANOS_INSTRUMENT (sys.getInstrumentation()->raiseOpenPtPEvent ( NANOS_WD_DOMAIN, (nanos_event_id_t) worker.getId(), 0, 0 ); )
   NANOS_INSTRUMENT (InstrumentationContextData *icd = worker.getInstrumentationContextData() );
   NANOS_INSTRUMENT (icd->setStartingWD(true) );
   //std::cout<<" associating thread\n"; 
   BaseThread &thread = createThread( worker );
   //std::cout<<" associated thread"<<&thread<<std::endl; 
   thread.associate();

   _threads.push_back( &thread );

   if ( !untieMain ) {
      worker.tieTo(thread);
   }

	//davidp
	///---int numevents=3;
        ///---int Events[numevents];// = {PAPI_TOT_INS, PAPI_TOT_CYC};
        ///---Events[0]= PAPI_TOT_INS;
        ///---Events[1]= PAPI_TOT_CYC;
	//Events[2]= PAPI_FDV_INS;
	///---Events[2]= PAPI_FP_OPS;

        ///---int retval=0;
	//sys.eventsets= new int[sys.getNumThreads()];
	sys.eventsets.reserve(sys.getNumThreads());
	
	//for (int i=0; i<sys.getNumThreads(); i++){
		//std::cout<<"two\n";
		sys.eventsets[getId()]=PAPI_NULL;
		myThread->eventse=PAPI_NULL;
	//}

   	//std::cout<<" Getting papi\n"; 
	//PAPI_shutdown();
	///---if (not PAPI_is_initialized()) {
	
                ///---if((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT ){
                ///---        std::cout<<"PAPI ERROR init (master): " << retval << "\n";
                ///---}else{ //std::cout<<"PAPI ERROR no error init (master): " << retval << "\n"; 
		///---	 }
                ///---	if((retval = PAPI_thread_init(pthread_self)) != PAPI_OK){
                ///---        std::cout<<"PAPI ERROR thread init (master): " << retval << "\n";
                ///---	}else{
			//std::cout<<"PAPI ERROR no errorno error thread init (master): " << retval << "\n";
		///---	}
	
        ///---}else{
              ///---  std::cout<<"PAPI ERROR already initialized (master): \n";
        ///---}
	//for(int i=0; i< sys.getNumThreads(); i++){

///////////

		std::cout << "one thread m" << "\n";
		///---if ((retval = PAPI_register_thread() ) != PAPI_OK){
                ///---        std::cout<<"PAPI ERROR register (thread): " << retval << "\n";
                ///---}else{//std::cout<<"PAPI ERROR no error register (thread): "  <<getId()<< "\n"; 
		///--- }

        	///---if ((retval=PAPI_create_eventset(&(myThread->eventse))) != PAPI_OK){
                ///---	std::cout<<"PAPI ERROR create set (master): " << retval << "\n";
        	///---}else{std::cout<<"PAPI ERROR no error create (master): "  <<getId()<< "\n";  }
        	///---if ((retval=PAPI_add_events(myThread->eventse, (int*)Events, numevents)) != PAPI_OK){
                ///---	std::cout<<"PAPI ERROR adding event to set (master!): " << retval << "\n";
        	///---}else{std::cout<<"PAPI ERROR no error add (master): "  <<getId()<< "\n";  }
        	///---if((retval=PAPI_start(myThread->eventse )) != PAPI_OK){
                ///---	std::cout<<"PAPI ERROR starting set (master): " << retval << "\n";
        	///---}else{std::cout<<"PAPI ERROR no error start (master): "  <<getId()<< "\n";  }

		//Events[0]= PAPI_L1_TCM;
        	//Events[1]= PAPI_L2_TCM;
        	//Events[2]= PAPI_L3_TCM;
	//}

	//\davidp

	/*if((retval=PAPI_start(sys.eventsets[0] )) != PAPI_OK){
                        std::cout<<"PAPI ERROR starting set (master): " << retval <<"0" << "\n";
                }else{std::cout<<"PAPI ERROR no error start (master): "  <<"0"<< "\n";  }
	
	if((retval=PAPI_start(sys.eventsets[1] )) != PAPI_OK){
                        std::cout<<"PAPI ERROR starting set (master): " << retval << "\n";
                }else{std::cout<<"PAPI ERROR no error start (master): "  <<"1"<< "\n";  }
*/

   return thread;
}

void ProcessingElement::stopAll ()
{
   ThreadList::iterator it;
   BaseThread *thread;

   for ( it = _threads.begin(); it != _threads.end(); it++ ) {
      thread = *it;
      if ( thread->getId() == 0) continue; /* Protection for master thread */
      thread->wakeup();
      thread->stop();
      thread->signal();
      thread->join();
      if ( thread->hasTeam() )
         thread->leaveTeam();
   }
}

void* ProcessingElement::getAddress( WorkDescriptor &wd, uint64_t tag, nanos_sharing_t sharing )
{
   void *actualTag = (void *) ( sharing == NANOS_PRIVATE ? (char *)wd.getData() + (unsigned long)tag : (void *)tag );
   return actualTag;
}

void ProcessingElement::copyTo( WorkDescriptor& wd, void *dst, uint64_t tag, nanos_sharing_t sharing, size_t size )
{
   void *actualTag = (void *) ( sharing == NANOS_PRIVATE ? (char *)wd.getData() + (unsigned long)tag : (void *)tag );
   // FIXME: should this be done by using the local copeir of the device?
   memcpy( dst, actualTag, size );
}

BaseThread* ProcessingElement::getFirstRunningThread()
{
   ThreadList::iterator it;
   for ( it = _threads.begin(); it != _threads.end(); it++ ) {
      if ( (*it)->hasTeam() && !(*it)->isTaggedToSleep() )
         return (*it);
   }
   return NULL;
}

BaseThread* ProcessingElement::getFirstStoppedThread()
{
   ThreadList::iterator it;
   for ( it = _threads.begin(); it != _threads.end(); it++ ) {
      if ( !(*it)->hasTeam() || (*it)->isTaggedToSleep() )
         return (*it);
   }
   return NULL;
}
