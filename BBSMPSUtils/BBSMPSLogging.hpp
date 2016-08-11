/*########################################################################
Copyright (c) 2014-2016, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.

Created by Geoffrey Oxberry (oxberry1@llnl.gov, goxberry@gmail.com),
Lluis-Miquel Munguia Conejero (lluis.munguia@gatech.edu), and Deepak
Rajan (rajan3@llnl.gov). LLNL-CODE-699387. All rights reserved.

This file is part of PIPS-SBB. For details, see
https://github.com/llnl/PIPS-SBB.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (as
published by the Free Software Foundation) version 2.1, February 1999.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
conditions of the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
########################################################################*/
// ----------------------------------------------------------------------------
/**
   File: BBSMPSLogging.hpp

   Description: Logging class responsible for centralizing the output of the
   				BB tree and activities.

   Limitations: Tags are a bit limited at the moment, but it will be expanded
   				in the future.

*/
// ----------------------------------------------------------------------------

#ifndef BBSMPSLOGGING_H
#define BBSMPSLOGGING_H
#include "PIPSLogging.hpp"



class BBSMPSLogging: public PIPSLogging
{
private:

  BBSMPSLogging() {};
public:

  //logging functions, however, MACROS defined at the end of the file are faster
  //and easier to use
  static void BBSMPSAppLogSeverity(severity_level lvl, std::string msg)
  {

    BOOST_LOG_SEV(PIPSLogging::g_sev_log, lvl) << logging::add_value("BBTREE", "on")<< msg;
  }
  static void BBSMPSAlgLogSeverity(severity_level lvl, std::string msg)
  {
    BOOST_LOG_SEV(PIPSLogging::g_sev_log, lvl) << logging::add_value("BBTREE", "on")<<logging::add_value("ALGORITHM", "on") << msg;
  }

  static void init_logging(int level)
  {
    using boost::shared_ptr;
    //create a text output sink:
    typedef sinks::synchronous_sink< sinks::text_ostream_backend > text_sink;
    shared_ptr< text_sink > pSink(new text_sink);

    // Here synchronous_sink is a sink frontend that performs thread synchronization
    // before passing log records to the backend; this makes backend easier to implement

    text_sink::locked_backend_ptr pBackend = pSink->locked_backend();

    // Add the stream corresponding to the console
    shared_ptr< std::ostream > pStream(&std::clog, boost::null_deleter());
    pBackend->add_stream(pStream);

      // Add the sink to the logging library
    logging::core::get()->add_sink(pSink);
    //set the formatter for the sink's output
    pSink->set_formatter(expr::stream << expr::if_(expr::has_attr("BBTREE"))
			 [
			  expr::stream << "BBSMPSTree: "
			   ]
			 <<
			 expr::if_(expr::has_attr("ALGORITHM"))
			 [
			  expr::stream << "  " << expr::smessage//expr::attr< std::string >("ALGORITHM")
			   ]
			 .else_
			 [
			  expr::stream
			  << "["
			  << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%H:%M:%S.%f")
			  << "] [" << expr::attr< severity_level >("Severity") << "] "
			  << expr::smessage
			  ]);

    pSink->set_filter(expr::attr< severity_level >("Severity").or_default(info) >= level);

    //add a time stamp
    attrs::local_clock TimeStamp;
    logging::core::get()->add_global_attribute("TimeStamp", TimeStamp);
  };
};

#define BBSMPS_APP_LOG_SEV(lvl) BOOST_LOG_SEV(PIPSLogging::g_sev_log, lvl)<< logging::add_value("BBTREE", "on")
#define BBSMPS_ALG_LOG_SEV(lvl) BOOST_LOG_SEV(PIPSLogging::g_sev_log, lvl) << logging::add_value("BBTREE", "on")<< logging::add_value("ALGORITHM", "on")


#endif
