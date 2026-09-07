#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void finalize_distances_parallel(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *precomputedArgsAndUpperTri, __global uchar *finalDistances, __private int N)
{
  long l_22, l_21, l_51, l_50, l_80, l_79, l_38, l_37, l_67, l_66, l_97, l_96, l_109, l_108, l_10, l_9; 
  ulong ul_28, ul_27, ul_26, ul_25, ul_32, ul_31, ul_29, ul_35, ul_33, ul_40, ul_12, ul_16, ul_14, ul_13, ul_17, ul_23, ul_0, ul_4, ul_3, ul_1, ul_112, ul_110, ul_116, ul_115, ul_114, ul_113, ul_119, ul_117, ul_91, ul_90, ul_89, ul_93, ul_100, ul_99, ul_104, ul_103, ul_101, ul_74, ul_73, ul_84, ul_83, ul_81, ul_87, ul_86, ul_85, ul_60, ul_58, ul_57, ul_64, ul_62, ul_61, ul_71, ul_70, ul_69, ul_44, ul_42, ul_41, ul_45, ul_52, ul_56, ul_55, ul_54; 
  int i_7, i_8, i_77, i_78, i_75, i_76, i_18, i_19, i_20, i_88, i_30, i_94, i_95, i_36, i_105, i_106, i_46, i_107, i_49, i_47, i_48, i_120, i_59, i_65, i_2, i_5, i_6; 
  bool b_39, b_102, b_72, b_11, b_43, b_15, b_98, b_68; 
  double d_92, d_24, d_118, d_53, d_34, d_82, d_63, d_111; 

  // BLOCK 0
  ul_0  =  (ulong) precomputedArgsAndUpperTri;
  ul_1  =  (ulong) finalDistances;
  i_2  =  get_global_size(0);
  ul_3  =  ul_0 + 32L;
  ul_4  =  ul_1 + 32L;
  i_5  =  get_global_id(0);
  // BLOCK 1 MERGES [0 22 ]
  i_6  =  i_5;
  for(;i_6 < 4;)
  {
    // BLOCK 2
    i_7  =  i_6 << 2;
    i_8  =  i_7 + 3;
    l_9  =  (long) i_8;
    l_10  =  l_9 << 3;
    b_11  =  i_6 == 0;
    if(b_11)
    {
      // BLOCK 3
      ul_12  =  *((__global ulong *) ul_4);
      ul_13  =  ul_1 + ul_12;
      ul_14  =  ul_13 + l_10;
      *((__global double *) ul_14)  =  0.0;
    }  // B3
    else
    {
      // BLOCK 4
      b_15  =  i_6 < 0;
      if(b_15)
      {
        // BLOCK 5
        ul_16  =  *((__global ulong *) ul_3);
        ul_17  =  ul_0 + ul_16;
        i_18  =  i_6 << 4;
        i_19  =  i_18 + i_7;
        i_20  =  i_19 + 7;
        l_21  =  (long) i_20;
        l_22  =  l_21 << 3;
        ul_23  =  ul_17 + l_22;
        d_24  =  *((__global double *) ul_23);
        ul_25  =  *((__global ulong *) ul_4);
        ul_26  =  ul_1 + ul_25;
        ul_27  =  ul_26 + l_10;
        *((__global double *) ul_27)  =  d_24;
      }  // B5
      else
      {
        // BLOCK 6
        ul_28  =  *((__global ulong *) ul_4);
        ul_29  =  ul_1 + ul_28;
        i_30  =  i_6 + 3;
        ul_31  =  i_30 & 4294967295UL;
        ul_32  =  ul_31 << 3;
        ul_33  =  ul_29 + ul_32;
        d_34  =  *((__global double *) ul_33);
        ul_35  =  ul_29 + l_10;
        *((__global double *) ul_35)  =  d_34;
      }  // B6
    }  // B4
    // BLOCK 7 MERGES [3 6 5 ]
    i_36  =  i_7 + 4;
    l_37  =  (long) i_36;
    l_38  =  l_37 << 3;
    b_39  =  i_6 == 1;
    if(b_39)
    {
      // BLOCK 8
      ul_40  =  *((__global ulong *) ul_4);
      ul_41  =  ul_1 + ul_40;
      ul_42  =  ul_41 + l_38;
      *((__global double *) ul_42)  =  0.0;
    }  // B8
    else
    {
      // BLOCK 9
      b_43  =  i_6 < 1;
      if(b_43)
      {
        // BLOCK 10
        ul_44  =  *((__global ulong *) ul_3);
        ul_45  =  ul_0 + ul_44;
        i_46  =  i_7 + 1;
        i_47  =  i_46 << 2;
        i_48  =  i_47 + i_46;
        i_49  =  i_48 + 7;
        l_50  =  (long) i_49;
        l_51  =  l_50 << 3;
        ul_52  =  ul_45 + l_51;
        d_53  =  *((__global double *) ul_52);
        ul_54  =  *((__global ulong *) ul_4);
        ul_55  =  ul_1 + ul_54;
        ul_56  =  ul_55 + l_38;
        *((__global double *) ul_56)  =  d_53;
      }  // B10
      else
      {
        // BLOCK 11
        ul_57  =  *((__global ulong *) ul_4);
        ul_58  =  ul_1 + ul_57;
        i_59  =  i_6 + 7;
        ul_60  =  i_59 & 4294967295UL;
        ul_61  =  ul_60 << 3;
        ul_62  =  ul_58 + ul_61;
        d_63  =  *((__global double *) ul_62);
        ul_64  =  ul_58 + l_38;
        *((__global double *) ul_64)  =  d_63;
      }  // B11
    }  // B9
    // BLOCK 12 MERGES [8 11 10 ]
    i_65  =  i_7 + 5;
    l_66  =  (long) i_65;
    l_67  =  l_66 << 3;
    b_68  =  i_6 == 2;
    if(b_68)
    {
      // BLOCK 13
      ul_69  =  *((__global ulong *) ul_4);
      ul_70  =  ul_1 + ul_69;
      ul_71  =  ul_70 + l_67;
      *((__global double *) ul_71)  =  0.0;
    }  // B13
    else
    {
      // BLOCK 14
      b_72  =  i_6 < 2;
      if(b_72)
      {
        // BLOCK 15
        ul_73  =  *((__global ulong *) ul_3);
        ul_74  =  ul_0 + ul_73;
        i_75  =  i_7 + 2;
        i_76  =  i_75 << 2;
        i_77  =  i_76 + i_75;
        i_78  =  i_77 + 7;
        l_79  =  (long) i_78;
        l_80  =  l_79 << 3;
        ul_81  =  ul_74 + l_80;
        d_82  =  *((__global double *) ul_81);
        ul_83  =  *((__global ulong *) ul_4);
        ul_84  =  ul_1 + ul_83;
        ul_85  =  ul_84 + l_67;
        *((__global double *) ul_85)  =  d_82;
      }  // B15
      else
      {
        // BLOCK 16
        ul_86  =  *((__global ulong *) ul_4);
        ul_87  =  ul_1 + ul_86;
        i_88  =  i_6 + 11;
        ul_89  =  i_88 & 4294967295UL;
        ul_90  =  ul_89 << 3;
        ul_91  =  ul_87 + ul_90;
        d_92  =  *((__global double *) ul_91);
        ul_93  =  ul_87 + l_67;
        *((__global double *) ul_93)  =  d_92;
      }  // B16
    }  // B14
    // BLOCK 17 MERGES [13 16 15 ]
    i_94  =  i_2 + i_6;
    i_95  =  i_7 + 6;
    l_96  =  (long) i_95;
    l_97  =  l_96 << 3;
    b_98  =  i_6 == 3;
    if(b_98)
    {
      // BLOCK 18
      ul_99  =  *((__global ulong *) ul_4);
      ul_100  =  ul_1 + ul_99;
      ul_101  =  ul_100 + l_97;
      *((__global double *) ul_101)  =  0.0;
    }  // B18
    else
    {
      // BLOCK 19
      b_102  =  i_6 < 3;
      if(b_102)
      {
        // BLOCK 20
        ul_103  =  *((__global ulong *) ul_3);
        ul_104  =  ul_0 + ul_103;
        i_105  =  i_8 << 2;
        i_106  =  i_105 + i_8;
        i_107  =  i_106 + 7;
        l_108  =  (long) i_107;
        l_109  =  l_108 << 3;
        ul_110  =  ul_104 + l_109;
        d_111  =  *((__global double *) ul_110);
        ul_112  =  *((__global ulong *) ul_4);
        ul_113  =  ul_1 + ul_112;
        ul_114  =  ul_113 + l_97;
        *((__global double *) ul_114)  =  d_111;
      }  // B20
      else
      {
        // BLOCK 21
        ul_115  =  *((__global ulong *) ul_4);
        ul_116  =  ul_1 + ul_115;
        ul_117  =  ul_116 + 144L;
        d_118  =  *((__global double *) ul_117);
        ul_119  =  ul_116 + l_97;
        *((__global double *) ul_119)  =  d_118;
      }  // B21
    }  // B19
    // BLOCK 22 MERGES [18 20 21 ]
    i_120  =  i_94;
    i_6  =  i_120;
  }  // B22
  // BLOCK 23
  return;
}  //  kernel
