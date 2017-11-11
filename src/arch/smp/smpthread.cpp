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
#include "os.hpp"
#include "smpprocessor.hpp"
#include "schedule.hpp"
#include "debug.hpp"
#include "system.hpp"
#include <iostream>
#include <sched.h>
#include <unistd.h>
#include "smp_ult.hpp"
#include "instrumentation.hpp"

#include <fstream>
#include <sstream> //davidp
//#include "system_decl.hpp"

#include <sys/time.h>

using namespace nanos;
using namespace nanos::ext;

pthread_mutex_t SMPThread::_mutexWait = PTHREAD_MUTEX_INITIALIZER;



template <typename T>
std::string to_string(T value)
{
	std::ostringstream os ;
	os << value ;
	return os.str() ;
}

void * smp_bootthread ( void *arg )
{
//---std::cout << "smp-bootthread ini" <<std::endl;
   SMPThread *self = static_cast<SMPThread *>( arg );
	
	////davidp
	///---int numevents=3;
	///---int Events[numevents];// = {PAPI_TOT_INS, PAPI_TOT_CYC};
        //long long values[NUM_EVENTS];
	///---Events[0]= PAPI_TOT_INS;
        ///---Events[1]= PAPI_TOT_CYC;
        //Events[2]= PAPI_FDV_INS;
        ///---Events[2]= PAPI_FP_OPS;
	//int EventSet=PAPI_NULL;
        ///---int retval=0;
        /*if (not PAPI_is_initialized()) {
                if((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT ){
                        std::cout<<"PAPI ERROR init bootthread.smpthread: " << retval << "\n";
                }else{ std::cout<<"PAPI ERROR no error init boothread.smpthread: " << retval << "\n";  }
                if((retval = PAPI_thread_init(pthread_self)) != PAPI_OK){
                        std::cout<<"PAPI ERROR thread init bootthread.smpthread: " << retval << "\n";
                }
        }else{
                std::cout<<"PAPI ERROR already initialized bootthread.smpthread: \n";
        }*/
	self->eventse=PAPI_NULL;
//---std::cout << "smp-bootthread before papi" <<std::endl;
	///---if ((retval=PAPI_create_eventset(&self->eventse)) != PAPI_OK){
        ///---        std::cout<<"PAPI ERROR create set (thread): " << retval << "\n";
        ///---}else{ //std::cout<<"PAPI ERROR no error starting set (thread): " << retval << "\n"; 
	///---}
        ///---if ((retval=PAPI_add_events(self->eventse, Events, numevents)) != PAPI_OK){
        ///---        std::cout<<"PAPI ERROR adding event to set (thread): " << retval << "\n";
        ///---}else{ //std::cout<<"PAPI ERROR no error starting set (thread): " << retval << "\n";  
	///---}
        ///---if((retval=PAPI_start(self->eventse)) != PAPI_OK){
        ///---        std::cout<<"PAPI ERROR starting set (thread): " << retval << "\n";
        ///---}else{// std::cout<<"PAPI ERROR no error starting set (thread): " << retval << "\n"; 
	///---}
	//////davidp

//---std::cout << "smp-bootthread after papi" <<std::endl;
   self->run();
//---std::cout << "smp-bootthread after papi run" <<std::endl;

   NANOS_INSTRUMENT ( static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary(); )
//---std::cout << "smp-bootthread after nstrumentation" <<std::endl;
   NANOS_INSTRUMENT ( static nanos_event_key_t cpuid_key = ID->getEventKey("cpuid"); )
   NANOS_INSTRUMENT ( nanos_event_value_t cpuid_value =  (nanos_event_value_t) 0; )
   NANOS_INSTRUMENT ( sys.getInstrumentation()->raisePointEvents(1, &cpuid_key, &cpuid_value); )

	//davidp
	/*if ((retval=PAPI_stop(EventSet, values)) != PAPI_OK){
                std::cout<<"PAPI ERROR stoping set: " << retval << "\n";
        }*/
	//davidp
/*	
	Events[0]= PAPI_TOT_INS;
        Events[1]= PAPI_TOT_CYC;
        Events[2]= PAPI_FDV_INS;
	sys.eventsets[myThread->getId()]=PAPI_NULL;

                std::cout << "one thread" << "\n";
                if ((retval=PAPI_create_eventset(&sys.eventsets[myThread->getId()])) != PAPI_OK){
                        std::cout<<"PAPI ERROR create set (thread): " << retval << "\n";
                }else{std::cout<<"PAPI ERROR no error create (thread): "  <<myThread->getId()<< "\n";  }
                if ((retval=PAPI_add_events(sys.eventsets[myThread->getId()], Events, numevents)) != PAPI_OK){
                        std::cout<<"PAPI ERROR adding event to set (thread): " << retval << "\n";
                }else{std::cout<<"PAPI ERROR no error add (thread): "  <<myThread->getId()<< "\n";  }
                if((retval=PAPI_start(sys.eventsets[myThread->getId()] )) != PAPI_OK){
                        std::cout<<"PAPI ERROR starting set (thread): " << retval << "\n";
                }else{std::cout<<"PAPI ERROR no error start (thread): "  <<myThread->getId()<< "\n";  }

*/
//---std::cout << "smp-bootthread end" <<std::endl;
   pthread_exit ( 0 );
   // We should never get here!
   return NULL;
}

// TODO: detect at configure
#ifndef PTHREAD_STACK_MIN
#define PTHREAD_STACK_MIN 16384
#endif

void SMPThread::start ()
{

   pthread_attr_t attr;
   pthread_attr_init(&attr);

   // user-defined stack size
   if ( _stackSize > 0 ) {
     if ( _stackSize < PTHREAD_STACK_MIN ) {
       warning("specified thread stack too small, adjusting it to minimum size");
       _stackSize = PTHREAD_STACK_MIN;
     }

     if (pthread_attr_setstacksize( &attr, _stackSize ) )
       warning("couldn't set pthread stack size stack");
   }

   if ( pthread_create( &_pth, &attr, smp_bootthread, this ) )
      fatal( "couldn't create thread" );

   if ( pthread_cond_init( &_condWait, NULL ) < 0 )
      fatal( "couldn't create pthread condition wait" );

	//int numevents=3;
        //int Events[numevents];	
	//int retval=0;
	/*Events[0]= PAPI_TOT_INS;
        Events[1]= PAPI_TOT_CYC;
        Events[2]= PAPI_FDV_INS;
        sys.eventsets[getId()]=PAPI_NULL;

	eventse=PAPI_NULL;*/
/////----->	std::cout<<"starting  and papizing"<<this<<std::endl;
/*
                std::cout << "one thread" << "\n";
		if ((retval = PAPI_register_thread() ) != PAPI_OK){
			std::cout<<"PAPI ERROR register (thread): " << retval << "\n";
		}else{std::cout<<"PAPI ERROR no error register (thread): "  <<getId()<< "\n";  }	

                if ((retval=PAPI_create_eventset(&eventse)) != PAPI_OK){
                        std::cout<<"PAPI ERROR create set (thread): " << retval << "\n";
                }else{std::cout<<"PAPI ERROR no error create (thread): "  <<getId()<< "\n";  }
                if ((retval=PAPI_add_events(eventse, Events, numevents)) != PAPI_OK){
                        std::cout<<"PAPI ERROR adding event to set (thread): " << retval << "\n";
                }else{std::cout<<"PAPI ERROR no error add (thread): "  <<getId()<< "\n";  }
                if((retval=PAPI_start(eventse )) != PAPI_OK){
                        std::cout<<"PAPI ERROR starting set (thread): " << retval << "\n";
                }else{std::cout<<"PAPI ERROR no error start (thread): "  <<getId()<< "\n";  }
*/	
}

void SMPThread::runDependent ()
{
   WD &work = getThreadWD();
   setCurrentWD( work );

   SMPDD &dd = ( SMPDD & ) work.activateDevice( SMP );

   dd.getWorkFct()( work.getData() );
}

void SMPThread::join ()
{
   if ( pthread_cond_destroy( &_condWait ) < 0 )
      fatal( "couldn't destroy pthread condition wait" );

   pthread_join( _pth,NULL );
   joined();
}

void SMPThread::bind( void )
{
   int cpu_id = getCpuId();

   cpu_set_t cpu_set;
   CPU_ZERO( &cpu_set );
   CPU_SET( cpu_id, &cpu_set );
   verbose( " Binding thread " << getId() << " to cpu " << cpu_id );
   OS::bindThread( &cpu_set );

   NANOS_INSTRUMENT ( static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary(); )
   NANOS_INSTRUMENT ( static nanos_event_key_t cpuid_key = ID->getEventKey("cpuid"); )
   NANOS_INSTRUMENT ( nanos_event_value_t cpuid_value =  (nanos_event_value_t) getCpuId() + 1; )
   NANOS_INSTRUMENT ( sys.getInstrumentation()->raisePointEvents(1, &cpuid_key, &cpuid_value); )
}

void SMPThread::yield()
{
   if (sched_yield() != 0)
      warning("sched_yield call returned an error");
}

void SMPThread::wait()
{
   NANOS_INSTRUMENT ( static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary(); )
   NANOS_INSTRUMENT ( static nanos_event_key_t cpuid_key = ID->getEventKey("cpuid"); )
   NANOS_INSTRUMENT ( nanos_event_value_t cpuid_value = (nanos_event_value_t) 0; )
   NANOS_INSTRUMENT ( sys.getInstrumentation()->raisePointEvents(1, &cpuid_key, &cpuid_value); )

   pthread_mutex_lock( &_mutexWait );

   if ( isTaggedToSleep() ) {
      ThreadTeam *team = getTeam();

      if ( hasNextWD() ) {
         WD *next = getNextWD();
         next->untie();
         team->getSchedulePolicy().queue( this, *next );
      }
      fatal_cond( hasNextWD(), "Can't sleep a thread with more than 1 WD in its local queue" );

      if ( team != NULL ) {
         team->removeThread( getTeamId() );
         leaveTeam();
      }
      pthread_cond_wait( &_condWait, &_mutexWait );
   }

   pthread_mutex_unlock( &_mutexWait );

   NANOS_INSTRUMENT ( if ( sys.getBinding() ) { cpuid_value = (nanos_event_value_t) getCpuId() + 1; } )
   NANOS_INSTRUMENT ( if ( !sys.getBinding() && sys.isCpuidEventEnabled() ) { cpuid_value = (nanos_event_value_t) sched_getcpu() + 1; } )
   NANOS_INSTRUMENT ( sys.getInstrumentation()->raisePointEvents(1, &cpuid_key, &cpuid_value); )
}

void SMPThread::signal()
{
   pthread_cond_signal( &_condWait );
}

void SMPThread::sleep()
{
   pthread_mutex_lock( &_mutexWait );
   BaseThread::sleep();
   pthread_mutex_unlock( &_mutexWait );
}

void SMPThread::wakeup()
{
   pthread_mutex_lock( &_mutexWait );
   BaseThread::wakeup();
   pthread_mutex_unlock( &_mutexWait );
}

void SMPThread::changePrefetcher(int reg, int pnum){
/*struct timeval t0, t1;
gettimeofday(&t0, NULL);*/
//FILE *in;
//char buff[512];
///adapt prefetcher register configuration number to vendor's format
std::string regc;
if(reg==0)//disabled 001
	regc="0x01";
else if(reg==1)//default
	regc="0x00";
else if(reg==2)//shallowest
	regc="0x02";
else if(reg==3)//deepest
	regc="0x07";
else if(reg==4)//prefetch on stores
	regc="0x08";
else if(reg==5)//aggressive
	regc="0x1f";
//std::cout<<"changed prefetch: "<<reg<< std::endl;
//else if(reg==4)//stride-N
//	reg=16;
//mapping for POWER7
//pnum=pnum*4;
pnum=((((pnum % 32)*4) + ((pnum % 32)/8)) % 32);
//std::cout << "number proc is: " << pnum << std::endl;
//-std::cout << "changing core: "<<pnum<<" value " <<reg <<std::endl;
///pnum=0;
//"echo "+to_string(reg)+" > /sys/devices/system/cpu/cpu"+spnum+"/dscr
std::stringstream ss;
//std::stringstream ss1;
ss << pnum;
//ss1 << reg;
std::fstream f;
////std::string file="/sys/devices/system/cpu/cpu"+ss.str()+"/dscr";
////const char * c = file.c_str();
char c[50];
sprintf(c, "/sys/devices/system/cpu/cpu%d/dscr", pnum);
//char file2[] = "/sys/devices/system/cpu/cpu0/dscr";
f.open( c , std::fstream::in | std::fstream::out | std::fstream::trunc ); //1
f << regc;//ss1.str();																										//2	
f.close();																																//3
/*gettimeofday(&t1, NULL); double time0=t0.tv_sec+(t0.tv_usec/1000000.0);
double time1=t1.tv_sec+(t1.tv_usec/1000000.0); printf("Seconds %f \n", time1-time0);*/
checkPrefetcher(pnum);
}

std::string SMPThread::checkPrefetcher(int pnum){
/*struct timeval t0, t1;
gettimeofday(&t0, NULL);*/

///FILE *in;
///char buff[50];
//mapping power7
//pnum=pnum*4
pnum=((((pnum % 32)*4) + ((pnum % 32)/8)) % 32);
//pnum=0;
	//std::string spnum=to_string(pnum);
	//std::string one="cat /sys/devices/system/cpu/cpu";
	//std::string two="/dscr";
	///std::stringstream ss;
	///ss << pnum;
	char tw[50];
	//sprintf(tw, "cat /sys/devices/system/cpu/cpu%d/dscr", pnum);
	sprintf(tw, "/sys/devices/system/cpu/cpu%d/dscr", pnum);
	///std::string file="/sys/devices/system/cpu/cpu"+ss.str()+"/dscr";
	///const char * tw = file.c_str();
	//std::string towrite =one.append(spnum);
	//towrite =towrite.append(two);
	//char * tw=strdup(towrite.c_str());
	//if(!(in = popen("cat /sys/devices/system/cpu/cpu"+spnum+"/dscr", "r"))){
	std::ifstream myReadFile;
  myReadFile.open(tw);
	char output[100];
 	if (myReadFile.is_open()) {
 	while (!myReadFile.eof()) {
    myReadFile >> output;
    }
	}
	myReadFile.close();

	//std::cout<< "dscr value: "<< std::string(output) <<std::endl;
	/*gettimeofday(&t1, NULL); double time0=t0.tv_sec+(t0.tv_usec/1000000.0);
double time1=t1.tv_sec+(t1.tv_usec/1000000.0); printf("Seconds %f \n", time1-time0);*/
	//return std::stoi(output);
	/////return atoi(output);///when returning an int
	return std::string(output);
/*
	if(!(in = popen(tw, "r"))){ //std::cout << "ERROR popen"<< std::endl ;//return 1; 
	}
  while(fgets(buff, sizeof(buff), in)!=NULL){ //std::cout << "(form smpthread) checking fgets value of the DSCR is: " << buff << std::endl;
  }
  pclose(in);
*/
/*
	if(buff!=NULL)
		return atoi(buff);
	else 
		return 0;
*/
}

void SMPThread::idle()
{

}

// This is executed in between switching stacks
void SMPThread::switchHelperDependent ( WD *oldWD, WD *newWD, void *oldState  )
{
   SMPDD & dd = ( SMPDD & )oldWD->getActiveDevice();
   dd.setState( (intptr_t *) oldState );
}

bool SMPThread::inlineWorkDependent ( WD &wd, WD &oldwd )
{//more than one thread
	
 	///---int ret=0;
	///---if ((ret=PAPI_read( eventse , values)) != PAPI_OK){
	///---	std::cout<<"PAPI ERROR inlineWorkDependent read oldwd: ret"<< ret <<" Id" <<getId() <<"\n";
	///---}else{//std::cout<<"PAPI ERROR no error inlineWorkDependent read oldwd: "<< " oldid"<< oldwd.getId()<<" ret" <<ret <<" Id"  <<getId()<<"\n"; 
	///---}
   	oldwd.addvalues(values);

   // Now the WD will be inminently run
   wd.start(WD::IsNotAUserLevelThread);
   //std::cout<<"Running task inline" << getId() <<std::endl;
   SMPDD &dd = ( SMPDD & )wd.getActiveDevice();

   NANOS_INSTRUMENT ( static nanos_event_key_t key = sys.getInstrumentation()->getInstrumentationDictionary()->getEventKey("user-code") );
   NANOS_INSTRUMENT ( nanos_event_value_t val = wd.getId() );
   NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateAndBurst ( NANOS_RUNNING, key, val ) );

   ///---if (PAPI_reset( eventse ) != PAPI_OK){ std::cout<<"PAPI ERROR reseting Id"<< getId() <<"\n"; }

   ( dd.getWorkFct() )( wd.getData() ); //wd2 call

   ///---if ((ret=PAPI_read( eventse , values)) != PAPI_OK){
   ///---   std::cout<<"PAPI ERROR inlineWorkDependent read oldwd: "<< ret <<" "<<getId()<< "\n";
   ///---}else{//std::cout<<"PAPI ERROR no error inlineWorkDependent read wd: "<< " wdId"<< wd.getId()<<" ret" <<ret<<" Id" << getId() << "\n";
	///--- }
   wd.addvalues(values);

   ///---if (PAPI_reset( eventse ) != PAPI_OK){ std::cout<<"PAPI ERROR reseting wdId" <<wd.getId()<<" Id" << getId()<< "\n"; }

   NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateAndBurst ( key, val ) );

   return true;
}

void SMPThread::switchTo ( WD *wd, SchedulerHelper *helper, WD &oldwd )
{
	///---int ret=0;
	///---if ((ret=PAPI_read( eventse ,values)) != PAPI_OK){
        ///---        std::cout<<"PAPI ERROR switchTo read oldwd "<< ret << "\n";
        ///---}else{// std::cout<<"PAPI ERROR no error switchTo read oldwd "<< ret <<" " << oldwd.getId()<<" " << values[0] << " " << values[1] <<" "<< getId() <<"\n"; 
///--- }
        oldwd.addvalues(values);	

	///---if (PAPI_reset( eventse  != PAPI_OK)) { std::cout<<"PAPI ERROR reseting "  << getId()<<"\n"; }

	//if ((ret=PAPI_start(sys.eventsets )) != PAPI_OK){std::cout<<"PAPI ERROR start! "<<ret<<"\n";}
   // wd MUST have an active SMP Deevice when it gets here
   //std::cout<<"Running task switchTo\n"<< getId()<<std::endl;
   ensure( wd->hasActiveDevice(),"WD has no active SMP device" );
   SMPDD &dd = ( SMPDD & )wd->getActiveDevice();
   ensure( dd.hasStack(), "DD has no stack for ULT");

   ::switchStacks(
       ( void * ) getCurrentWD(),
       ( void * ) wd,
       ( void * ) dd.getState(),
       ( void * ) helper );

	///---if (PAPI_reset( eventse ) != PAPI_OK){ std::cout<<"PAPI ERROR reseting\n"; }
	//std::cout<<"Running task switchTo end\n"<<std::endl;

}

void SMPThread::exitTo ( WD *wd, SchedulerHelper *helper, WD &oldwd)
{

   //davidp: code for reading PAPI counters is in wordescriptor finish()
   // wd MUST have an active SMP Device when it gets here
   ensure( wd->hasActiveDevice(),"WD has no active SMP device" );
   SMPDD &dd = ( SMPDD & )wd->getActiveDevice();
   ensure( dd.hasStack(), "DD has no stack for ULT");

   //std::cout<<"Running task exitTo "<< getId() <<std::endl;
   ///---if (PAPI_reset( eventse ) != PAPI_OK){ std::cout<<"PAPI ERROR reseting "<< getId() <<"\n"; }//reset before starting new wd
   //TODO: optimize... we don't really need to save a context in this case
   ::switchStacks(
      ( void * ) getCurrentWD(),
      ( void * ) wd,
      ( void * ) dd.getState(),
      ( void * ) helper );
	//non reachable code
   	std::cout<<"Running task exitTo end\n"<<std::endl;
	
}

