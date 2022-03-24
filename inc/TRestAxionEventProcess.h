/*************************************************************************
 * This file is part of the REST software framework.                     *
 *                                                                       *
 * Copyright (C) 2016 GIFNA/TREX (University of Zaragoza)                *
 * For more information see http://gifna.unizar.es/trex                  *
 *                                                                       *
 * REST is free software: you can redistribute it and/or modify          *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * REST is distributed in the hope that it will be useful,               *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have a copy of the GNU General Public License along with   *
 * REST in $REST_PATH/LICENSE.                                           *
 * If not, see http://www.gnu.org/licenses/.                             *
 * For the list of contributors see $REST_PATH/CREDITS.                  *
 *************************************************************************/

#ifndef RestCore_TRestAxionEventProcess
#define RestCore_TRestAxionEventProcess

#include "TRestEventProcess.h"

/// A base class for any REST event process
class TRestAxionEventProcess : public TRestEventProcess {
   private:
    // unused datamember, set as private
    ///< not used, keep for compatibility
    // TRestEvent* fInputEvent = nullptr;  //!
    ///< not used, keep for compatibility
    // TRestEvent* fOutputEvent = nullptr;  //!

   protected:
    // utils
    void BeginPrintProcess();
    void EndPrintProcess();

   public:
    // process running methods
    /// To be executed at the beginning of the run (outside event loop)
    virtual void InitProcess() {}
    /// Begin of event process, preparation work. Called right before ProcessEvent()
    void BeginOfEventProcess(TRestEvent* evInput = nullptr);
    /// End of event process. Nothing to do. Called directly after ProcessEvent()
    void EndOfEventProcess(TRestEvent* evInput = nullptr);

    TRestAxionEventProcess();
    ~TRestAxionEventProcess();

    ClassDef(TRestAxionEventProcess, 1);
};
#endif
