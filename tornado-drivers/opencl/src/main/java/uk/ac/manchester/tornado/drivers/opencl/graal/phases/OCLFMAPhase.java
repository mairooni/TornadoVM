/*
 * Copyright (c) 2020, 2025, APT Group, Department of Computer Science,
 * The University of Manchester. All rights reserved.
 * Copyright (c) 2009, 2017, Oracle and/or its affiliates. All rights reserved.
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 *
 * This code is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 only, as
 * published by the Free Software Foundation.
 *
 * This code is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * version 2 for more details (a copy is included in the LICENSE file that
 * accompanied this code).
 *
 * You should have received a copy of the GNU General Public License version
 * 2 along with this work; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */
package uk.ac.manchester.tornado.drivers.opencl.graal.phases;

import java.util.Optional;

import org.graalvm.compiler.nodes.GraphState;
import org.graalvm.compiler.nodes.StructuredGraph;
import org.graalvm.compiler.nodes.ValueNode;
import org.graalvm.compiler.nodes.calc.AddNode;
import org.graalvm.compiler.nodes.calc.MulNode;
import org.graalvm.compiler.phases.Phase;

import jdk.vm.ci.meta.JavaKind;
import uk.ac.manchester.tornado.api.TornadoDeviceContext;
import uk.ac.manchester.tornado.drivers.opencl.graal.nodes.AddHalfNode;
import uk.ac.manchester.tornado.drivers.opencl.graal.nodes.MultHalfNode;
import uk.ac.manchester.tornado.drivers.opencl.graal.nodes.OCLFMANode;

public class OCLFMAPhase extends Phase {

    private TornadoDeviceContext deviceContext;

    public OCLFMAPhase(TornadoDeviceContext deviceContext) {
        this.deviceContext = deviceContext;
    }

    /**
     * Instrinsics in OpenCL:
     * https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/fma.html
     * 
     * 
     * From the OpenCL documentation:
     * 
     * The generic type name gentype is used to indicate that the function can take
     * float, float2, float4, float8, or float16 as the type for the arguments. For
     * any specific use of these function, the actual type has to be the same for
     * all arguments and the return type.
     *
     * If extended with cl_khr_fp64, generic type name gentype may indicate double
     * and double{2|4|8|16} as arguments and return values. If extended with
     * cl_khr_fp16, generic type name gentype may indicate half and half{2|4|8|16}
     * as arguments and return values.
     * 
     * @param x
     *            valueNode
     * @return returns true if ValueNode is either Float of Double
     */
    private boolean isValidType(ValueNode x) {
        return (x.getStackKind() == JavaKind.Float || x.getStackKind() == JavaKind.Double);
    }

    @Override
    public Optional<NotApplicable> notApplicableTo(GraphState graphState) {
        return ALWAYS_APPLICABLE;
    }

    private void applyFMA(ValueNode addNode, ValueNode mulNode, ValueNode x, ValueNode y, StructuredGraph graph) {
        ValueNode z = (ValueNode) addNode.inputs().filter(node -> !node.equals(mulNode)).first();
        OCLFMANode oclFMA = new OCLFMANode(x, y, z);
        graph.addOrUnique(oclFMA);
        mulNode.removeUsage(addNode);
        if (mulNode.hasNoUsages()) {
            mulNode.safeDelete();
        }
        addNode.replaceAtUsages(oclFMA);
        addNode.safeDelete();
    }

    @Override
    protected void run(StructuredGraph graph) {
        boolean isNVIDIAGPU = deviceContext.getDeviceName().contains("NVIDIA");
        // FMA is not currently available for NVIDIA GPUs in OpenCL
        if (!isNVIDIAGPU) {
            graph.getNodes().filter(AddHalfNode.class).forEach(addHalfNode -> {
                MultHalfNode multHalfNode = null;
                if (addHalfNode.getX() instanceof MultHalfNode) {
                    multHalfNode = (MultHalfNode) addHalfNode.getX();
                } else if (addHalfNode.getY() instanceof MultHalfNode) {
                    multHalfNode = (MultHalfNode) addHalfNode.getY();
                }
                if (multHalfNode != null) {
                    ValueNode x = multHalfNode.getX();
                    ValueNode y = multHalfNode.getY();
                    applyFMA(addHalfNode, multHalfNode, x, y, graph);
                }
            });
        }

        graph.getNodes().filter(AddNode.class).forEach(addNode -> {
            MulNode mulNode = null;
            if (addNode.getX() instanceof MulNode) {
                mulNode = (MulNode) addNode.getX();
            } else if (addNode.getY() instanceof MulNode) {
                mulNode = (MulNode) addNode.getY();
            }
            if (mulNode != null) {
                ValueNode x = mulNode.getX();
                ValueNode y = mulNode.getY();
                if (isValidType(x) && isValidType(y)) {
                    applyFMA(addNode, mulNode, x, y, graph);
                }
            }
        });

    }
}
