//MODIFIED 2 WIDE VLIW

///////////////////////////////////////////
// wallypipelinedcore.sv
//
// Written: David_Harris@hmc.edu 9 January 2021
// Modified:
//
// Purpose: Pipelined RISC-V Processor
//
// Documentation: RISC-V System on Chip Design
//
// A component of the CORE-V-WALLY configurable RISC-V project.
// https://github.com/openhwgroup/cvw
//
// Copyright (C) 2021-23 Harvey Mudd College & Oklahoma State University
//
// SPDX-License-Identifier: Apache-2.0 WITH SHL-2.1
//
// Licensed under the Solderpad Hardware License v 2.1 (the “License”); you may not use this file
// except in compliance with the License, or, at your option, the Apache License version 2.0. You
// may obtain a copy of the License at
//
// https://solderpad.org/licenses/SHL-2.1/
//
// Unless required by applicable law or agreed to in writing, any work distributed under the
// License is distributed on an “AS IS” BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
// either express or implied. See the License for the specific language governing permissions
// and limitations under the License.
////////////////////////////////////////////////////////////////////////////////////////////////

module wallypipelinedcore import cvw::*; #(parameter cvw_t P) (
   input  logic                  clk, reset,
   // Privileged
   input  logic                  MTimerInt, MExtInt, SExtInt, MSwInt,
   input  logic [63:0]           MTIME_CLINT,
   // Bus Interface
   input  logic [P.AHBW-1:0]     HRDATA,
   input  logic                  HREADY, HRESP,
   output logic                  HCLK, HRESETn,
   output logic [P.PA_BITS-1:0]  HADDR,
   output logic [P.AHBW-1:0]     HWDATA,
   output logic [P.XLEN/8-1:0]   HWSTRB,
   output logic                  HWRITE,
   output logic [2:0]            HSIZE,
   output logic [2:0]            HBURST,
   output logic [3:0]            HPROT,
   output logic [1:0]            HTRANS,
   output logic                  HMASTLOCK,
   input  logic                  ExternalStall
);

  logic                          StallF, StallD, StallE, StallM, StallW;
  logic                          FlushD, FlushE, FlushM, FlushW;
  logic                          TrapM, RetM;

  //  signals that must connect through DP
  logic                          IntDivE, W64E;
  logic                          CSRReadM, CSRWriteM, PrivilegedM;
  logic [1:0]                    AtomicM;
  logic [P.XLEN-1:0]             ForwardedSrcAE, ForwardedSrcBE;
  logic [P.XLEN-1:0]             SrcAM;
  logic [2:0]                    Funct3E;
  logic [31:0]                   InstrD;
  logic [31:0]                   InstrM, InstrOrigM;
  logic [P.XLEN-1:0]             PCSpillF, PCE, PCLinkE;
  logic [P.XLEN-1:0]             PCM, PCSpillM;
  logic [P.XLEN-1:0]             CSRReadValW, MDUResultW;
  logic [P.XLEN-1:0]             EPCM, TrapVectorM;
  logic [1:0]                    MemRWE;
  logic [1:0]                    MemRWM;
  logic                          InstrValidD, InstrValidE, InstrValidM;
  logic                          InstrMisalignedFaultM;
  logic                          IllegalBaseInstrD, IllegalFPUInstrD, IllegalIEUFPUInstrD;
  logic                          InstrPageFaultF, LoadPageFaultM, StoreAmoPageFaultM;
  logic                          LoadMisalignedFaultM, LoadAccessFaultM;
  logic                          StoreAmoMisalignedFaultM, StoreAmoAccessFaultM;
  logic                          InvalidateICacheM, FlushDCacheM;
  logic                          PCSrcE;
  logic                          CSRWriteFenceM;
  logic                          DivBusyE;
  logic                          StructuralStallD;
  logic                          LoadStallD;
  logic                          StoreStallD;
  logic                          SquashSCW;
  logic                          MDUActiveE;                      // Mul/Div instruction being executed
  logic                          ENVCFG_ADUE;                     // HPTW A/D Update enable
  logic                          ENVCFG_PBMTE;                    // Page-based memory type enable
  logic [3:0]                    ENVCFG_CBE;                      // Cache Block operation enables
  logic [3:0]                    CMOpM;                           // 1: cbo.inval; 2: cbo.flush; 4: cbo.clean; 8: cbo.zero
  logic                          IFUPrefetchE, LSUPrefetchM;      // instruction / data prefetch hints


  // floating point unit signals
  logic [2:0]                    FRM_REGW;
  logic [4:0]                    RdE, RdM, RdW;
  logic                          FPUStallD;
  logic                          FWriteIntE;
  logic [P.FLEN-1:0]             FWriteDataM;
  logic [P.XLEN-1:0]             FIntResM;
  logic [P.XLEN-1:0]             FCvtIntResW;
  logic                          FCvtIntW;
  logic                          FDivBusyE;
  logic                          FRegWriteM;
  logic                          FpLoadStoreM;
  logic [4:0]                    SetFflagsM;
  logic [P.XLEN-1:0]             FIntDivResultW;

  // memory management unit signals
  logic                          ITLBWriteF;
  logic                          ITLBMissOrUpdateAF;
  logic [P.XLEN-1:0]             SATP_REGW;
  logic                          STATUS_MXR, STATUS_SUM, STATUS_MPRV;
  logic [1:0]                    STATUS_MPP, STATUS_FS;
  logic [1:0]                    PrivilegeModeW;
  logic [P.XLEN-1:0]             PTE;
  logic [1:0]                    PageType;
  logic                          sfencevmaM;
  logic                          SelHPTW;

  // PMA checker signals
  /* verilator lint_off UNDRIVEN */ // these signals are undriven in configurations without a privileged unit
  var logic [P.PA_BITS-3:0]      PMPADDR_ARRAY_REGW[P.PMP_ENTRIES-1:0];
  var logic [7:0]                PMPCFG_ARRAY_REGW[P.PMP_ENTRIES-1:0];
  /* verilator lint_on UNDRIVEN */

  // IMem stalls
  logic                          IFUStallF;
  logic                          LSUStallM;

  // cpu lsu interface
  logic [2:0]                    Funct3M;
  logic [P.XLEN-1:0]             IEUAdrE;
  logic [P.XLEN-1:0]             WriteDataM;
  logic [P.XLEN-1:0]             IEUAdrM;
  logic [P.XLEN-1:0]             IEUAdrxTvalM;
  logic [P.LLEN-1:0]             ReadDataW;
  logic                          CommittedM;

  // AHB ifu interface
  logic [P.PA_BITS-1:0]          IFUHADDR;
  logic [2:0]                    IFUHBURST;
  logic [1:0]                    IFUHTRANS;
  logic [2:0]                    IFUHSIZE;
  logic                          IFUHWRITE;
  logic                          IFUHREADY;

  // AHB LSU interface
  logic [P.PA_BITS-1:0]          LSUHADDR;
  logic [P.XLEN-1:0]             LSUHWDATA;
  logic [P.XLEN/8-1:0]           LSUHWSTRB;
  logic                          LSUHWRITE;
  logic                          LSUHREADY;

  logic                          BPWrongE, BPWrongM;
  logic                          BPDirWrongM;
  logic                          BTAWrongM;
  logic                          RASPredPCWrongM;
  logic                          IClassWrongM;
  logic [3:0]                    IClassM;
  logic                          InstrAccessFaultF, HPTWInstrAccessFaultF, HPTWInstrPageFaultF;
  logic [2:0]                    LSUHSIZE;
  logic [2:0]                    LSUHBURST;
  logic [1:0]                    LSUHTRANS;

  logic                          DCacheMiss;
  logic                          DCacheAccess;
  logic                          ICacheMiss;
  logic                          ICacheAccess;
  logic                          BigEndianM;
  logic                          FCvtIntE;
  logic                          CommittedF;
  logic                          BranchD, BranchE, JumpD, JumpE;
  logic                          DCacheStallM, ICacheStallF;
  logic                          wfiM, IntPendingM;

  //VLIW
  logic [31:0]                   VLIWInstr0D, VLIWInstr1D, VLIWInstr2D, VLIWInstr3D;
  logic [3:0]                    VLIWValidD;
  logic                          VLIWModeD;

  //VLIW Auxillary Signals
  logic [P.XLEN-1:0]             FCvtIntResW_0;
  logic [P.XLEN-1:0]             FIntDivResultW_0;
  logic [P.XLEN-1:0]             FIntResM_0;
  logic [P.XLEN-1:0]             ForwardedSrcAE_0, ForwardedSrcBE_0;
  // logic [P.XLEN-1:0]             IEUAdrE_0;
  logic [P.XLEN-1:0]             MDUResultW_0;
  logic [P.XLEN-1:0]             SrcAM_0;
  logic [P.XLEN-1:0]             WriteDataM_0;
  
  logic [P.FLEN-1:0]             FWriteDataM_0;
  logic [P.LLEN-1:0]             ReadDataW_0;

  logic [31:0]                   InstrM_0;

  logic [4:0]                    RdE_0, RdM_0, RdW_0;

  logic [3:0]                    CMOpM_0;

  logic [2:0]                    Funct3E_0;
  logic [2:0]                    Funct3M_0;
  
  logic [1:0]                    MemRWE_0, MemRWM_0;       
  logic [1:0]                    AtomicM_0;

  logic                          BranchD_0, BranchE_0, JumpD_0, JumpE_0;
  logic                          CSRWriteFenceM_0;
  logic                          FCvtIntE_0;
  logic                          FCvtIntW_0;
  logic                          FlushDCacheM_0;
  logic                          FpLoadStoreM_0;
  logic                          FWriteIntE_0;
  logic                          IFUPrefetchE_0, LSUPrefetchM_0;
  logic                          IllegalBaseInstrD_0;
  logic                          InstrValidM_0, InstrValidE_0, InstrValidD_0;
  logic                          IntDivE_0;
  logic                          MDUActiveE_0;
  logic                          SquashSCW_0;
  logic                          StructuralStallD_0;
  logic                          W64E_0;


  logic [P.XLEN-1:0]             FCvtIntResW_1;
  logic [P.XLEN-1:0]             FIntDivResultW_1;
  logic [P.XLEN-1:0]             FIntResM_1;
  logic [P.XLEN-1:0]             ForwardedSrcAE_1, ForwardedSrcBE_1;
  // logic [P.XLEN-1:0]             IEUAdrE_1;
  logic [P.XLEN-1:0]             MDUResultW_1;
  logic [P.XLEN-1:0]             SrcAM_1;
  logic [P.XLEN-1:0]             WriteDataM_1;
  
  logic [P.FLEN-1:0]             FWriteDataM_1;
  logic [P.LLEN-1:0]             ReadDataW_1;

  logic [31:0]                   InstrM_1;

  logic [4:0]                    RdE_1, RdM_1, RdW_1;

  logic [3:0]                    CMOpM_1;

  logic [2:0]                    Funct3E_1;
  logic [2:0]                    Funct3M_1;
  
  logic [1:0]                    MemRWE_1, MemRWM_1;       
  logic [1:0]                    AtomicM_1;

  logic                          BranchD_1, BranchE_1, JumpD_1, JumpE_1;
  logic                          CSRWriteFenceM_1;
  logic                          FCvtIntE_1;
  logic                          FCvtIntW_1;
  logic                          FlushDCacheM_1;
  logic                          FpLoadStoreM_1;
  logic                          FWriteIntE_1;
  logic                          IFUPrefetchE_1, LSUPrefetchM_1;
  logic                          IllegalBaseInstrD_1;
  logic                          InstrValidM_1, InstrValidE_1, InstrValidD_1;
  logic                          IntDivE_1;
  logic                          MDUActiveE_1;
  logic                          SquashSCW_1;
  logic                          StructuralStallD_1;
  logic                          W64E_1;

  

  // instruction fetch unit: PC, branch prediction, instruction cache
  ifu #(P) ifu(.clk, .reset,
    .StallF, .StallD, .StallE, .StallM, .StallW, .FlushD, .FlushE, .FlushM, .FlushW,
    .InstrValidE, .InstrValidD,
    .BranchD, .BranchE, .JumpD, .JumpE, .ICacheStallF,
    // Fetch
    .HRDATA, .PCSpillF, .IFUHADDR,
    .IFUStallF, .IFUHBURST, .IFUHTRANS, .IFUHSIZE, .IFUHREADY, .IFUHWRITE,
    .ICacheAccess, .ICacheMiss,
    // Execute
    .PCLinkE, .PCSrcE, .IEUAdrE, .IEUAdrM, .PCE, .BPWrongE,  .BPWrongM,
    // Mem
    .CommittedF, .EPCM, .TrapVectorM, .RetM, .TrapM, .InvalidateICacheM, .CSRWriteFenceM,
    .InstrD, .InstrM, .InstrOrigM, .PCM, .PCSpillM, .IClassM, .BPDirWrongM,
    .BTAWrongM, .RASPredPCWrongM, .IClassWrongM,
    // Faults out
    .IllegalBaseInstrD, .IllegalFPUInstrD, .InstrPageFaultF, .IllegalIEUFPUInstrD, .InstrMisalignedFaultM,
    // mmu management
    .PrivilegeModeW, .PTE, .PageType, .SATP_REGW, .STATUS_MXR, .STATUS_SUM, .STATUS_MPRV,
    .STATUS_MPP, .ENVCFG_PBMTE, .ENVCFG_ADUE, .ITLBWriteF, .sfencevmaM, .ITLBMissOrUpdateAF,
    // pmp/pma (inside mmu) signals.
    .PMPCFG_ARRAY_REGW,  .PMPADDR_ARRAY_REGW, .InstrAccessFaultF,
    // VLIW mode
    .VLIWInstr0D,  .VLIWInstr1D,  .VLIWInstr2D,  .VLIWInstr3D,  .VLIWValidD,  .VLIWModeD);

  // integer execution unit: integer register file, datapath and controller
  ieu #(P) 
    ieu_0(.clk, .reset,
      // Decode Stage interface
      .InstrD, .STATUS_FS, .ENVCFG_CBE, .IllegalIEUFPUInstrD, .IllegalBaseInstrD,
      // Execute Stage interface
      .PCE, .PCLinkE, .FWriteIntE(FWriteIntE_0), .FCvtIntE(FCvtIntE_0), .IEUAdrE, .IntDivE(IntDivE_0), .W64E(W64E_0),
      .Funct3E(Funct3E_0), .ForwardedSrcAE(ForwardedSrcAE_0), .ForwardedSrcBE(ForwardedSrcBE_0), .MDUActiveE(MDUActiveE_0), .CMOpM(CMOpM_0), .IFUPrefetchE(IFUPrefetchE_0), .LSUPrefetchM(LSUPrefetchM_0),
      // Memory stage interface
      .SquashSCW(SquashSCW_0),  // from LSU
      .MemRWE(MemRWE_0),     // read/write control goes to LSU
      .MemRWM(MemRWM_0),     // read/write control goes to LSU
      .AtomicM(AtomicM_0),    // atomic control goes to LSU
      .WriteDataM(WriteDataM_0), // Write data to LSU
      .Funct3M(Funct3M_0),    // size and signedness to LSU
      .SrcAM(SrcAM_0),      // to privilege and fpu
      .RdE(RdE_0), .RdM(RdM_0), .FIntResM(FIntResM_0), .FlushDCacheM(FlushDCacheM_0),
      .BranchD(BranchD_0), .BranchE(BranchE_0), .JumpD(JumpD_0), .JumpE(JumpE_0),
      // Writeback stage
      .CSRReadValW, .MDUResultW(MDUResultW_0), .FIntDivResultW(FIntDivResultW_0), .RdW(RdW_0), .ReadDataW(ReadDataW_0[P.XLEN-1:0]),
      .InstrValidM(InstrValidM_0), .InstrValidE(InstrValidE_0), .InstrValidD(InstrValidD_0), .FCvtIntResW(FCvtIntResW_0), .FCvtIntW(FCvtIntW_0),
      // hazards
      .StallD, .StallE, .StallM, .StallW, .FlushD, .FlushE, .FlushM, .FlushW,
      .StructuralStallD(StructuralStallD_0), .LoadStallD, .StoreStallD, .PCSrcE,
      .CSRReadM, .CSRWriteM, .PrivilegedM, .CSRWriteFenceM(CSRWriteFenceM_0), .InvalidateICacheM);
  ieu #(P)
    ieu_1(.clk, .reset,
      // Decode Stage interface
      .InstrD, .STATUS_FS, .ENVCFG_CBE, .IllegalIEUFPUInstrD, .IllegalBaseInstrD,
      // Execute Stage interface
      .PCE, .PCLinkE, .FWriteIntE(FWriteIntE_1), .FCvtIntE(FCvtIntE_1), .IEUAdrE, .IntDivE(IntDivE_1), .W64E(W64E_1),
      .Funct3E(Funct3E_1), .ForwardedSrcAE(ForwardedSrcAE_1), .ForwardedSrcBE(ForwardedSrcBE_1), .MDUActiveE(MDUActiveE_1), .CMOpM(CMOpM_1), .IFUPrefetchE(IFUPrefetchE_1), .LSUPrefetchM(LSUPrefetchM_1),
      // Memory stage interface
      .SquashSCW(SquashSCW_1),  // from LSU
      .MemRWE(MemRWE_1),     // read/write control goes to LSU
      .MemRWM(MemRWM_1),     // read/write control goes to LSU
      .AtomicM(AtomicM_1),    // atomic control goes to LSU
      .WriteDataM(WriteDataM_1), // Write data to LSU
      .Funct3M(Funct3M_1),    // size and signedness to LSU
      .SrcAM(SrcAM_1),      // to privilege and fpu
      .RdE(RdE_1), .RdM(RdM_1), .FIntResM(FIntResM_1), .FlushDCacheM(FlushDCacheM_1),
      .BranchD(BranchD_1), .BranchE(BranchE_1), .JumpD(JumpD_1), .JumpE(JumpE_1),
      // Writeback stage
      .CSRReadValW, .MDUResultW(MDUResultW_1), .FIntDivResultW(FIntDivResultW_1), .RdW(RdW_1), .ReadDataW(ReadDataW_1[P.XLEN-1:0]),
      .InstrValidM(InstrValidM_1), .InstrValidE(InstrValidE_1), .InstrValidD(InstrValidD_1), .FCvtIntResW(FCvtIntResW_1), .FCvtIntW(FCvtIntW_1),
      // hazards
      .StallD, .StallE, .StallM, .StallW, .FlushD, .FlushE, .FlushM, .FlushW,
      .StructuralStallD(StructuralStallD_1), .LoadStallD, .StoreStallD, .PCSrcE,
      .CSRReadM, .CSRWriteM, .PrivilegedM, .CSRWriteFenceM(CSRWriteFenceM_1), .InvalidateICacheM);

  lsu #(P) 
    lsu_0(
      .clk, .reset, .StallM, .FlushM, .StallW, .FlushW,
      // CPU interface
      .MemRWE(MemRWE_0), .MemRWM(MemRWM_0), .Funct3M(Funct3M_0), .Funct7M(InstrM_0[31:25]), .AtomicM(AtomicM_0),
      .CommittedM, .DCacheMiss, .DCacheAccess, .SquashSCW(SquashSCW_0),
      .FpLoadStoreM(FpLoadStoreM_0), .FWriteDataM(FWriteDataM_0), .IEUAdrE, .IEUAdrM, .WriteDataM(WriteDataM_0),
      .ReadDataW(ReadDataW_0), .FlushDCacheM(FlushDCacheM_0), .CMOpM(CMOpM_0), .LSUPrefetchM(LSUPrefetchM_0),
      // connected to ahb (all stay the same)
      .LSUHADDR,  .HRDATA, .LSUHWDATA, .LSUHWSTRB, .LSUHSIZE,
      .LSUHBURST, .LSUHTRANS, .LSUHWRITE, .LSUHREADY,
      // connect to csr or privilege and stay the same.
      .PrivilegeModeW, .BigEndianM, // connects to csr
      .PMPCFG_ARRAY_REGW,           // connects to csr
      .PMPADDR_ARRAY_REGW,          // connects to csr
      // hptw keep i/o
      .SATP_REGW,                   // from csr
      .STATUS_MXR,                  // from csr
      .STATUS_SUM,                  // from csr
      .STATUS_MPRV,                 // from csr
      .STATUS_MPP,                  // from csr
      .ENVCFG_PBMTE,                // from csr
      .ENVCFG_ADUE,                 // from csr
      .sfencevmaM,                  // connects to privilege
      .DCacheStallM,                // connects to privilege
      .IEUAdrxTvalM,                // connects to privilege
      .LoadPageFaultM,              // connects to privilege
      .StoreAmoPageFaultM,          // connects to privilege
      .LoadMisalignedFaultM,        // connects to privilege
      .LoadAccessFaultM,            // connects to privilege
      .HPTWInstrAccessFaultF,       // connects to privilege
      .HPTWInstrPageFaultF,         // connects to privilege
      .StoreAmoMisalignedFaultM,    // connects to privilege
      .StoreAmoAccessFaultM,        // connects to privilege
      .PCSpillF, .ITLBMissOrUpdateAF, .PTE, .PageType, .ITLBWriteF, .SelHPTW,
      .LSUStallM);
  lsu #(P)
    lsu_1(
      .clk, .reset, .StallM, .FlushM, .StallW, .FlushW,
      // CPU interface
      .MemRWE(MemRWE_1), .MemRWM(MemRWM_1), .Funct3M(Funct3M_1), .Funct7M(InstrM_1[31:25]), .AtomicM(AtomicM_1),
      .CommittedM, .DCacheMiss, .DCacheAccess, .SquashSCW(SquashSCW_1),
      .FpLoadStoreM(FpLoadStoreM_1), .FWriteDataM(FWriteDataM_1), .IEUAdrE, .IEUAdrM, .WriteDataM(WriteDataM_1),
      .ReadDataW(ReadDataW_1), .FlushDCacheM(FlushDCacheM_1), .CMOpM(CMOpM_1), .LSUPrefetchM(LSUPrefetchM_1),
      // connected to ahb (all stay the same)
      .LSUHADDR,  .HRDATA, .LSUHWDATA, .LSUHWSTRB, .LSUHSIZE,
      .LSUHBURST, .LSUHTRANS, .LSUHWRITE, .LSUHREADY,
      // connect to csr or privilege and stay the same.
      .PrivilegeModeW, .BigEndianM, // connects to csr
      .PMPCFG_ARRAY_REGW,           // connects to csr
      .PMPADDR_ARRAY_REGW,          // connects to csr
      // hptw keep i/o
      .SATP_REGW,                   // from csr
      .STATUS_MXR,                  // from csr
      .STATUS_SUM,                  // from csr
      .STATUS_MPRV,                 // from csr
      .STATUS_MPP,                  // from csr
      .ENVCFG_PBMTE,                // from csr
      .ENVCFG_ADUE,                 // from csr
      .sfencevmaM,                  // connects to privilege
      .DCacheStallM,                // connects to privilege
      .IEUAdrxTvalM,                // connects to privilege
      .LoadPageFaultM,              // connects to privilege
      .StoreAmoPageFaultM,          // connects to privilege
      .LoadMisalignedFaultM,        // connects to privilege
      .LoadAccessFaultM,            // connects to privilege
      .HPTWInstrAccessFaultF,       // connects to privilege
      .HPTWInstrPageFaultF,         // connects to privilege
      .StoreAmoMisalignedFaultM,    // connects to privilege
      .StoreAmoAccessFaultM,        // connects to privilege
      .PCSpillF, .ITLBMissOrUpdateAF, .PTE, .PageType, .ITLBWriteF, .SelHPTW,
      .LSUStallM);

  // global stall and flush control
  hazard hzu(
    .BPWrongE, .CSRWriteFenceM, .RetM, .TrapM,
    .StructuralStallD,
    .LSUStallM, .IFUStallF,
    .FPUStallD, .ExternalStall,
    .DivBusyE, .FDivBusyE,
    .wfiM, .IntPendingM,
    // Stall & flush outputs
    .StallF, .StallD, .StallE, .StallM, .StallW,
    .FlushD, .FlushE, .FlushM, .FlushW);

  if(P.BUS_SUPPORTED) begin : ebu
    ebu #(P) ebu(// IFU connections
      .clk, .reset,
      // IFU interface
      .IFUHADDR, .IFUHBURST, .IFUHTRANS, .IFUHREADY, .IFUHSIZE,
      // LSU interface
      .LSUHADDR, .LSUHWDATA, .LSUHWSTRB, .LSUHSIZE, .LSUHBURST,
      .LSUHTRANS, .LSUHWRITE, .LSUHREADY,
      // BUS interface
      .HREADY, .HRESP, .HCLK, .HRESETn,
      .HADDR, .HWDATA, .HWSTRB, .HWRITE, .HSIZE, .HBURST,
      .HPROT, .HTRANS, .HMASTLOCK);
  end else begin
    assign {IFUHREADY, LSUHREADY, HCLK, HRESETn, HADDR, HWDATA,
            HWSTRB, HWRITE, HSIZE, HBURST, HPROT, HTRANS, HMASTLOCK} = '0;
  end



  // privileged unit
  if (P.ZICSR_SUPPORTED) begin:priv
    privileged #(P) priv(
      .clk, .reset,
      .FlushD, .FlushE, .FlushM, .FlushW, .StallD, .StallE, .StallM, .StallW,
      .CSRReadM, .CSRWriteM, .SrcAM, .PCM, .PCSpillM,
      .InstrM, .InstrOrigM, .CSRReadValW, .EPCM, .TrapVectorM,
      .RetM, .TrapM, .sfencevmaM, .InvalidateICacheM, .DCacheStallM, .ICacheStallF,
      .InstrValidM, .CommittedM, .CommittedF,
      .FRegWriteM, .LoadStallD, .StoreStallD,
      .BPDirWrongM, .BTAWrongM, .BPWrongM,
      .RASPredPCWrongM, .IClassWrongM, .DivBusyE, .FDivBusyE,
      .IClassM, .DCacheMiss, .DCacheAccess, .ICacheMiss, .ICacheAccess, .PrivilegedM,
      .InstrPageFaultF, .LoadPageFaultM, .StoreAmoPageFaultM,
      .InstrMisalignedFaultM, .IllegalIEUFPUInstrD,
      .LoadMisalignedFaultM, .StoreAmoMisalignedFaultM,
      .MTimerInt, .MExtInt, .SExtInt, .MSwInt,
      .MTIME_CLINT, .IEUAdrxTvalM, .SetFflagsM,
      .InstrAccessFaultF, .HPTWInstrAccessFaultF, .HPTWInstrPageFaultF, .LoadAccessFaultM, .StoreAmoAccessFaultM, .SelHPTW,
      .PrivilegeModeW, .SATP_REGW,
      .STATUS_MXR, .STATUS_SUM, .STATUS_MPRV, .STATUS_MPP, .STATUS_FS,
      .PMPCFG_ARRAY_REGW, .PMPADDR_ARRAY_REGW,
      .FRM_REGW, .ENVCFG_CBE, .ENVCFG_PBMTE, .ENVCFG_ADUE, .wfiM, .IntPendingM, .BigEndianM);
  end else begin
    assign {CSRReadValW, PrivilegeModeW,
            SATP_REGW, STATUS_MXR, STATUS_SUM, STATUS_MPRV, STATUS_MPP, STATUS_FS, FRM_REGW,
            // PMPCFG_ARRAY_REGW, PMPADDR_ARRAY_REGW,
            ENVCFG_CBE, ENVCFG_PBMTE, ENVCFG_ADUE,
            EPCM, TrapVectorM, RetM, TrapM,
            sfencevmaM, BigEndianM, wfiM, IntPendingM} = '0;
  end

  // multiply/divide unit
  if (P.ZMMUL_SUPPORTED) begin:mdu
    mdu #(P) 
    mdu_0(.clk, .reset, .StallM, .StallW, .FlushE, .FlushM, .FlushW,
      .ForwardedSrcAE(ForwardedSrcAE_0), .ForwardedSrcBE(ForwardedSrcBE_0),
      .Funct3E(Funct3E_0), .Funct3M(Funct3M_0), .IntDivE(IntDivE_0), .W64E(W64E_0), .MDUActiveE(MDUActiveE_0),
      .MDUResultW(MDUResultW_0), .DivBusyE);
  end else begin // no M instructions supported
    assign MDUResultW_0 = '0;
    assign DivBusyE   = 1'b0;
  end
  if (P.ZMMUL_SUPPORTED) begin:mdu
    mdu #(P) 
    mdu_1(.clk, .reset, .StallM, .StallW, .FlushE, .FlushM, .FlushW,
      .ForwardedSrcAE(ForwardedSrcAE_1), .ForwardedSrcBE(ForwardedSrcBE_1),
      .Funct3E(Funct3E_1), .Funct3M(Funct3M_1), .IntDivE(IntDivE_1), .W64E(W64E_1), .MDUActiveE(MDUActiveE_1),
      .MDUResultW(MDUResultW_1), .DivBusyE);
  end else begin // no M instructions supported
    assign MDUResultW_1 = '0;
    assign DivBusyE   = 1'b0;
  end

  // floating point unit
  if (P.F_SUPPORTED) begin:fpu
    fpu #(P) 
    fpu_0(
      .clk, .reset,
      .FRM_REGW,                           // Rounding mode from CSR
      .InstrD,                             // instruction from IFU
      .ReadDataW(ReadDataW_0[P.FLEN-1:0]),   // Read data from memory
      .ForwardedSrcAE(ForwardedSrcAE_0),                     // Integer input being processed (from IEU)
      .StallE, .StallM, .StallW,           // stall signals from HZU
      .FlushE, .FlushM, .FlushW,           // flush signals from HZU
      .RdE(RdE_0), .RdM(RdM_0), .RdW(RdW_0),                    // which FP register to write to (from IEU)
      .STATUS_FS,                          // is floating-point enabled?
      .FRegWriteM,                         // FP register write enable
      .FpLoadStoreM,
      .ForwardedSrcBE(ForwardedSrcBE_0),                     // Integer input for intdiv
      .Funct3E(Funct3E_0), .Funct3M(Funct3M_0), .IntDivE(IntDivE_0), .W64E(W64E_0), // Integer flags and functions
      .FPUStallD,                          // Stall the decode stage
      .FWriteIntE(FWriteIntE_0), .FCvtIntE(FCvtIntE_0),              // integer register write enable, conversion operation
      .FWriteDataM(FWriteDataM_0),                        // Data to be written to memory
      .FIntResM(FIntResM_0),                           // data to be written to integer register
      .FCvtIntResW(FCvtIntResW_0),                        // fp -> int conversion result to be stored in int register
      .FCvtIntW(FCvtIntW_0),                           // fpu result selection
      .FDivBusyE,                          // Is the divide/sqrt unit busy (stall execute stage)
      .IllegalFPUInstrD,                   // Is the instruction an illegal fpu instruction
      .SetFflagsM,                         // FPU flags (to privileged unit)
      .FIntDivResultW(FIntDivResultW_0));
  end else begin                           // no F_SUPPORTED or D_SUPPORTED; tie outputs low
    assign {FPUStallD, FWriteIntE_0, FCvtIntE_0, FIntResM_0, FCvtIntW_0, FRegWriteM,
            IllegalFPUInstrD, SetFflagsM, FpLoadStoreM,
            FWriteDataM_0, FCvtIntResW_0, FIntDivResultW_0, FDivBusyE} = '0;
  end
  if (P.F_SUPPORTED) begin:fpu
    fpu #(P) 
    fpu_1(
      .clk, .reset,
      .FRM_REGW,                           // Rounding mode from CSR
      .InstrD,                             // instruction from IFU
      .ReadDataW(ReadDataW_1[P.FLEN-1:0]),   // Read data from memory
      .ForwardedSrcAE(ForwardedSrcAE_1),                     // Integer input being processed (from IEU)
      .StallE, .StallM, .StallW,           // stall signals from HZU
      .FlushE, .FlushM, .FlushW,           // flush signals from HZU
      .RdE(RdE_1), .RdM(RdM_1), .RdW(RdW_1),                    // which FP register to write to (from IEU)
      .STATUS_FS,                          // is floating-point enabled?
      .FRegWriteM,                         // FP register write enable
      .FpLoadStoreM,
      .ForwardedSrcBE(ForwardedSrcBE_1),                     // Integer input for intdiv
      .Funct3E(Funct3E_1), .Funct3M(Funct3M_1), .IntDivE(IntDivE_1), .W64E(W64E_1), // Integer flags and functions
      .FPUStallD,                          // Stall the decode stage
      .FWriteIntE(FWriteIntE_1), .FCvtIntE(FCvtIntE_1),              // integer register write enable, conversion operation
      .FWriteDataM(FWriteDataM_1),                        // Data to be written to memory
      .FIntResM(FIntResM_1),                           // data to be written to integer register
      .FCvtIntResW(FCvtIntResW_1),                        // fp -> int conversion result to be stored in int register
      .FCvtIntW(FCvtIntW_1),                           // fpu result selection
      .FDivBusyE,                          // Is the divide/sqrt unit busy (stall execute stage)
      .IllegalFPUInstrD,                   // Is the instruction an illegal fpu instruction
      .SetFflagsM,                         // FPU flags (to privileged unit)
      .FIntDivResultW(FIntDivResultW_1));
  end else begin                           // no F_SUPPORTED or D_SUPPORTED; tie outputs low
    assign {FPUStallD, FWriteIntE_1, FCvtIntE_1, FIntResM_1, FCvtIntW_1, FRegWriteM,
            IllegalFPUInstrD, SetFflagsM, FpLoadStoreM,
            FWriteDataM_1, FCvtIntResW_1, FIntDivResultW_1, FDivBusyE} = '0;
  end

endmodule
