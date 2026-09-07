#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
float uk_ac_manchester_tornado_examples_compute_MatrixVectorRowMajor_matrixVectorRowMajorOptimizedQ4(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, ulong context, int localSize, ulong x, ulong weightsQ4, ulong weightScales, int n)
{
  char ch_45, ch_41, ch_155, ch_37, ch_49; 
  bool b_182; 
  long l_158, l_27, l_28, l_157, l_86, l_87, l_153, l_82, l_146, l_47, l_43, l_39, l_162, l_35, l_163, l_62, l_61, l_56, l_57, l_51, l_52, l_81, l_145, l_76, l_77, l_71, l_72, l_66, l_67; 
  half half_148, half_30; 
  ulong ul_44, ul_73, ul_40, ul_83, ul_147, ul_48, ul_78, ul_58, ul_154, ul_88, ul_53, ul_36, ul_68, ul_164, ul_63, ul_159, ul_29; 
  int i_156, i_26, i_25, i_24, i_152, i_23, i_151, i_22, i_150, i_34, i_33, i_161, i_32, i_31, i_13, i_141, i_12, i_140, i_11, i_139, i_10, i_9, i_137, i_8, i_136, i_7, i_135, i_6, i_134, i_21, i_149, i_19, i_18, i_17, i_16, i_144, i_15, i_143, i_14, i_142, i_120, i_119, i_5, i_133, i_4, i_132, i_3, i_131, i_2, i_1, i_108, i_107, i_104, i_103, i_116, i_115, i_112, i_111, i_92, i_91, i_100, i_99, i_96, i_95, i_75, i_70, i_85, i_80, i_60, i_187, i_184, i_55, i_65, i_172, i_171, i_42, i_168, i_167, i_38, i_181, i_180, i_179, i_50, i_178, i_177, i_46; 
  float f_20, f_84, f_89, f_74, f_138, f_79, f_101, f_165, f_98, f_105, f_169, f_102, f_166, f_93, f_90, f_160, f_97, f_94, f_117, f_114, f_121, f_185, f_54, f_118, f_183, f_109, f_173, f_106, f_170, f_176, f_113, f_110, f_174, f_175, f_69, f_130, f_124, f_188, f_125, f_122, f_186, f_59, f_123, f_64, f_128, f_129, f_126, f_127; 

  // BLOCK 0
  __local float adf_0[localSize];
  i_1  =  localSize << 3;
  i_2  =  get_group_id(0);
  i_3  =  n * i_2;
  i_4  =  i_3 >> 31;
  i_5  =  i_4 + i_3;
  i_6  =  i_5 >> 1;
  i_7  =  i_6 + 19;
  i_8  =  i_6 + 18;
  i_9  =  i_6 + 17;
  i_10  =  i_6 + 16;
  i_11  =  n >> 31;
  i_12  =  i_11 >> 27;
  i_13  =  i_12 + n;
  i_14  =  i_13 >> 5;
  i_15  =  i_14 * i_2;
  i_16  =  i_15 + 8;
  i_17  =  get_local_id(0);
  i_18  =  i_17 << 3;
  i_19  =  n + -7;
  // BLOCK 1 MERGES [0 2 ]
  f_20  =  0.0F;
  i_21  =  i_18;
  for(;i_21 < i_19;)
  {
    // BLOCK 2
    i_22  =  i_21 >> 31;
    i_23  =  i_22 >> 27;
    i_24  =  i_23 + i_21;
    i_25  =  i_24 >> 5;
    i_26  =  i_25 + i_16;
    l_27  =  (long) i_26;
    l_28  =  l_27 << 1;
    ul_29  =  weightScales + l_28;
    half_30  =  *((__global half *) ul_29);
    i_31  =  i_21 >> 31;
    i_32  =  i_31 + i_21;
    i_33  =  i_32 >> 1;
    i_34  =  i_33 + i_10;
    l_35  =  (long) i_34;
    ul_36  =  weightsQ4 + l_35;
    ch_37  =  *((__global char *) ul_36);
    i_38  =  i_9 + i_33;
    l_39  =  (long) i_38;
    ul_40  =  weightsQ4 + l_39;
    ch_41  =  *((__global char *) ul_40);
    i_42  =  i_8 + i_33;
    l_43  =  (long) i_42;
    ul_44  =  weightsQ4 + l_43;
    ch_45  =  *((__global char *) ul_44);
    i_46  =  i_7 + i_33;
    l_47  =  (long) i_46;
    ul_48  =  weightsQ4 + l_47;
    ch_49  =  *((__global char *) ul_48);
    i_50  =  i_21 + 4;
    l_51  =  (long) i_50;
    l_52  =  l_51 << 2;
    ul_53  =  x + l_52;
    f_54  =  *((__global float *) ul_53);
    i_55  =  i_21 + 5;
    l_56  =  (long) i_55;
    l_57  =  l_56 << 2;
    ul_58  =  x + l_57;
    f_59  =  *((__global float *) ul_58);
    i_60  =  i_21 + 6;
    l_61  =  (long) i_60;
    l_62  =  l_61 << 2;
    ul_63  =  x + l_62;
    f_64  =  *((__global float *) ul_63);
    i_65  =  i_21 + 7;
    l_66  =  (long) i_65;
    l_67  =  l_66 << 2;
    ul_68  =  x + l_67;
    f_69  =  *((__global float *) ul_68);
    i_70  =  i_21 + 8;
    l_71  =  (long) i_70;
    l_72  =  l_71 << 2;
    ul_73  =  x + l_72;
    f_74  =  *((__global float *) ul_73);
    i_75  =  i_21 + 9;
    l_76  =  (long) i_75;
    l_77  =  l_76 << 2;
    ul_78  =  x + l_77;
    f_79  =  *((__global float *) ul_78);
    i_80  =  i_21 + 10;
    l_81  =  (long) i_80;
    l_82  =  l_81 << 2;
    ul_83  =  x + l_82;
    f_84  =  *((__global float *) ul_83);
    i_85  =  i_21 + 11;
    l_86  =  (long) i_85;
    l_87  =  l_86 << 2;
    ul_88  =  x + l_87;
    f_89  =  *((__global float *) ul_88);
    f_90  =  convert_float((float) half_30);
    i_91  =  (int) ch_49;
    i_92  =  i_91 >> 4;
    f_93  =  (float) i_92;
    f_94  =  f_90 * f_93;
    i_95  =  i_91 << 28;
    i_96  =  i_95 >> 28;
    f_97  =  (float) i_96;
    f_98  =  f_97 * f_90;
    i_99  =  (int) ch_45;
    i_100  =  i_99 >> 4;
    f_101  =  (float) i_100;
    f_102  =  f_101 * f_90;
    i_103  =  i_99 << 28;
    i_104  =  i_103 >> 28;
    f_105  =  (float) i_104;
    f_106  =  f_105 * f_90;
    i_107  =  (int) ch_41;
    i_108  =  i_107 >> 4;
    f_109  =  (float) i_108;
    f_110  =  f_109 * f_90;
    i_111  =  i_107 << 28;
    i_112  =  i_111 >> 28;
    f_113  =  (float) i_112;
    f_114  =  f_113 * f_90;
    i_115  =  (int) ch_37;
    i_116  =  i_115 >> 4;
    f_117  =  (float) i_116;
    f_118  =  f_117 * f_90;
    i_119  =  i_115 << 28;
    i_120  =  i_119 >> 28;
    f_121  =  (float) i_120;
    f_122  =  f_121 * f_90;
    f_123  =  fma(f_122, f_54, f_20);
    f_124  =  fma(f_118, f_59, f_123);
    f_125  =  fma(f_114, f_64, f_124);
    f_126  =  fma(f_110, f_69, f_125);
    f_127  =  fma(f_106, f_74, f_126);
    f_128  =  fma(f_102, f_79, f_127);
    f_129  =  fma(f_98, f_84, f_128);
    f_130  =  fma(f_94, f_89, f_129);
    i_131  =  i_1 + i_21;
    f_20  =  f_130;
    i_21  =  i_131;
  }  // B2
  // BLOCK 3
  i_132  =  localSize << 1;
  i_133  =  i_17 << 1;
  i_134  =  i_11 >> 29;
  i_135  =  i_134 + n;
  i_136  =  i_135 & -8;
  i_137  =  i_133 + i_136;
  // BLOCK 4 MERGES [3 5 ]
  f_138  =  f_20;
  i_139  =  i_137;
  for(;i_139 < n;)
  {
    // BLOCK 5
    i_140  =  i_139 >> 31;
    i_141  =  i_140 >> 27;
    i_142  =  i_141 + i_139;
    i_143  =  i_142 >> 5;
    i_144  =  i_143 + i_16;
    l_145  =  (long) i_144;
    l_146  =  l_145 << 1;
    ul_147  =  weightScales + l_146;
    half_148  =  *((__global half *) ul_147);
    i_149  =  i_139 >> 31;
    i_150  =  i_149 + i_139;
    i_151  =  i_150 >> 1;
    i_152  =  i_151 + i_10;
    l_153  =  (long) i_152;
    ul_154  =  weightsQ4 + l_153;
    ch_155  =  *((__global char *) ul_154);
    i_156  =  i_139 + 4;
    l_157  =  (long) i_156;
    l_158  =  l_157 << 2;
    ul_159  =  x + l_158;
    f_160  =  *((__global float *) ul_159);
    i_161  =  i_139 + 5;
    l_162  =  (long) i_161;
    l_163  =  l_162 << 2;
    ul_164  =  x + l_163;
    f_165  =  *((__global float *) ul_164);
    f_166  =  convert_float((float) half_148);
    i_167  =  (int) ch_155;
    i_168  =  i_167 >> 4;
    f_169  =  (float) i_168;
    f_170  =  f_166 * f_169;
    i_171  =  i_167 << 28;
    i_172  =  i_171 >> 28;
    f_173  =  (float) i_172;
    f_174  =  f_173 * f_166;
    f_175  =  fma(f_174, f_160, f_138);
    f_176  =  fma(f_170, f_165, f_175);
    i_177  =  i_132 + i_139;
    f_138  =  f_176;
    i_139  =  i_177;
  }  // B5
  // BLOCK 6
  adf_0[i_17]  =  f_138;
  barrier(CLK_LOCAL_MEM_FENCE);
  i_178  =  localSize >> 31;
  i_179  =  i_178 + localSize;
  i_180  =  i_179 >> 1;
  // BLOCK 7 MERGES [6 11 ]
  i_181  =  i_180;
  for(;i_181 >= 1;)
  {
    // BLOCK 8
    b_182  =  i_17 < i_181;
    if(b_182)
    {
      // BLOCK 9
      f_183  =  adf_0[i_17];
      i_184  =  i_181 + i_17;
      f_185  =  adf_0[i_184];
      f_186  =  f_183 + f_185;
      adf_0[i_17]  =  f_186;
    }  // B9
    else
    {
      // BLOCK 10
    }  // B10
    // BLOCK 11 MERGES [10 9 ]
    barrier(CLK_LOCAL_MEM_FENCE);
    i_187  =  i_181 >> 1;
    i_181  =  i_187;
  }  // B11
  // BLOCK 12
  f_188  =  adf_0[0];
  return f_188;
}  //  kernel

#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void matrixVectorGenericQ4(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *x, __global uchar *output, __global uchar *weightsQ4, __global uchar *weightScales, __private int dim1, __private int dim0, __private int localWorkGroupSize)
{
  bool b_9, b_6; 
  long l_11, l_12; 
  ulong ul_4, ul_3, ul_2, ul_1, ul_0, ul_13; 
  int i_10, i_7, i_5; 
  float f_8; 

  // BLOCK 0
  ul_0  =  (ulong) context;
  ul_1  =  (ulong) x;
  ul_2  =  (ulong) output;
  ul_3  =  (ulong) weightsQ4;
  ul_4  =  (ulong) weightScales;
  i_5  =  get_group_id(0);
  b_6  =  i_5 < 2048;
  if(b_6)
  {
    // BLOCK 1
    i_7  =  get_local_id(0);
    f_8  =  uk_ac_manchester_tornado_examples_compute_MatrixVectorRowMajor_matrixVectorRowMajorOptimizedQ4(_kernel_context, _constant_region, _local_region, _atomics, ul_0, 32, ul_1, ul_3, ul_4, 8192);
    b_9  =  i_7 == 0;
    if(b_9)
    {
      // BLOCK 2
      i_10  =  i_5 + 4;
      l_11  =  (long) i_10;
      l_12  =  l_11 << 2;
      ul_13  =  ul_2 + l_12;
      *((__global float *) ul_13)  =  f_8;
      return;
    }  // B2
    else
    {
      // BLOCK 3
      return;
    }  // B3
    // BLOCK 4 MERGES [3 2 ]
    return;
    else
    {
      // BLOCK 5
      return;
    }  // B5
  }  //  kernel
