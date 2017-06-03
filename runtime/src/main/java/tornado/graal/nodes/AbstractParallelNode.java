/* 
 * Copyright 2012 James Clarkson.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package tornado.graal.nodes;

import com.oracle.graal.compiler.common.type.StampFactory;
import com.oracle.graal.graph.NodeClass;
import com.oracle.graal.nodeinfo.NodeInfo;
import com.oracle.graal.nodes.ValueNode;
import com.oracle.graal.nodes.calc.FloatingNode;
import jdk.vm.ci.meta.JavaKind;

@NodeInfo
public abstract class AbstractParallelNode extends FloatingNode implements Comparable<AbstractParallelNode> {

    public static final NodeClass<AbstractParallelNode> TYPE = NodeClass
            .create(AbstractParallelNode.class);

    protected int index;
    @Input
    protected ValueNode value;

    protected AbstractParallelNode(NodeClass<? extends AbstractParallelNode> type, int index, ValueNode value) {
        super(type, StampFactory.forKind(JavaKind.Int));
        assert stamp != null;
        this.index = index;
        this.value = value;
    }

    public ValueNode value() {
        return value;
    }

    public int index() {
        return index;
    }

    @Override
    public int compareTo(AbstractParallelNode o) {
        return Integer.compare(index, o.index);
    }

}