/*
 * This file is part of Tornado: A heterogeneous programming framework:
 * https://github.com/beehive-lab/tornadovm
 *
 * Copyright (c) 2021, 2025, APT Group, Department of Computer Science,
 * School of Engineering, The University of Manchester. All rights reserved.
 * Copyright (c) 2009-2021, Oracle and/or its affiliates. All rights reserved.
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
package uk.ac.manchester.tornado.drivers.spirv.graal.phases;

import java.util.Optional;

import org.graalvm.compiler.nodes.GraphState;
import org.graalvm.compiler.nodes.StructuredGraph;
import org.graalvm.compiler.nodes.ValueNode;
import org.graalvm.compiler.nodes.calc.AddNode;
import org.graalvm.compiler.nodes.calc.MulNode;
import org.graalvm.compiler.phases.Phase;

import jdk.vm.ci.meta.JavaKind;
import uk.ac.manchester.tornado.drivers.spirv.graal.nodes.AddHalfNode;
import uk.ac.manchester.tornado.drivers.spirv.graal.nodes.MultHalfNode;
import uk.ac.manchester.tornado.drivers.spirv.graal.nodes.SPIRVFMANode;

public class SPIRVFMAPhase extends Phase {

    public Optional<NotApplicable> notApplicableTo(GraphState graphState) {
        return ALWAYS_APPLICABLE;
    }

    private boolean isValidType(ValueNode x) {
        return (x.getStackKind() == JavaKind.Float || x.getStackKind() == JavaKind.Double);
    }

    private void applyFMA(ValueNode addNode, ValueNode mulNode, ValueNode x, ValueNode y, StructuredGraph graph) {
        ValueNode z = (ValueNode) addNode.inputs().filter(node -> !node.equals(mulNode)).first();
        SPIRVFMANode oclFMA = new SPIRVFMANode(x, y, z);
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
