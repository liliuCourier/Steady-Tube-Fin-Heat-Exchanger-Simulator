function [heatPaths, pdropPaths,predecessors_in,predecessors_out] = buildPaths(A, N)
    [n_nodes, n_pipes] = size(A);
    
    % ---------- 解析每根管的流入/流出节点 ----------
    in_node  = zeros(1, n_pipes);
    out_node = zeros(1, n_pipes);
    inlet_pipes  = [];
    outlet_pipes = [];
    
    for p = 1:n_pipes
        col = A(:, p);
        in_idx  = find(col == -1);
        out_idx = find(col == 1);
        
        if ~isempty(out_idx) && ~isempty(in_idx)
            out_node(p) = out_idx;
            in_node(p)  = in_idx;
        elseif ~isempty(in_idx) && isempty(out_idx)
            in_node(p)  = in_idx;
            out_node(p) = 0;
            inlet_pipes(end+1) = p;
        elseif isempty(in_idx) && ~isempty(out_idx)
            out_node(p) = out_idx;
            in_node(p)  = 0;
            outlet_pipes(end+1) = p;
        else
            error('管道 %d 无连接', p);
        end
    end
    
    % ---------- 构建管间依赖关系 ----------
    % 对于每个节点，找出所有流入管 (in_pipes) 和流出管 (out_pipes)
    node_in_pipes  = cell(1, n_nodes);
    node_out_pipes = cell(1, n_nodes);
    for p = 1:n_pipes
        if in_node(p) > 0
            node_in_pipes{in_node(p)}(end+1) = p;
        end
        if out_node(p) > 0
            node_out_pipes{out_node(p)}(end+1) = p;
        end
    end
    
    % 构建有向无环图：管 -> 依赖的管（前驱）
    % 规则：对于每个节点，流出管依赖于该节点的所有流入管（因为汇合时需要所有上游信息）
    predecessors_in = cell(1, n_pipes);
    for node = 1:n_nodes
        ins  = node_in_pipes{node};
        outs = node_out_pipes{node};
        if isempty(ins) || isempty(outs)
            continue;
        end
        % 每个流出管依赖所有流入管
        for j = 1:length(outs)
            predecessors_in{outs(j)} = unique([predecessors_in{outs(j)}, ins]);
        end
    end

    predecessors_out = cell(1, n_pipes);
    for node = 1:n_nodes
        ins  = node_in_pipes{node};
        outs = node_out_pipes{node};
        if isempty(ins) || isempty(outs)
            continue;
        end
        % 每个流出管依赖所有流入管
        for j = 1:length(ins)
            predecessors_out{ins(j)} = unique([predecessors_out{ins(j)}, outs]);
        end
    end

    % ---------- 拓扑排序（标准Kahn算法） ----------
    % 广度优先算法本质上就是个遍历算法，这样来写
    indeg = zeros(1, n_pipes);
    for i = 1:n_pipes
        indeg(i) = length(predecessors_in{i});
    end
    % 队列初始为所有入度为0的管（即进口管）
    queue = find(indeg == 0);
    isIn = [];
    notIn = 1:n_pipes;
    isIn(1,end+1) = queue(1);
    notIn( notIn == queue(1)) = [];

    heatPaths = [];
    heatPaths(1,1) = queue(1);

    % 开始遍历
    while size(notIn,2) > 0
        for i = 1:size(notIn,2)
            tube = notIn(i);
            if all(ismember(predecessors_in{tube},isIn)) || ismember(tube,queue)
                isIn(1,end+1) = tube;
                heatPaths(1,end+1) = tube;
                notIn(notIn == tube) = [];
                break
            end
        end
    end

    % ---------- 生成压降路径（环路优先） ----------
    pdropPaths = [];
    road_num = size(N,2);

    for i = 1:road_num
        r1 =  find(N(:,i)==1)';
        r2 =  find(N(:,i)==-1)';

        % 排序r1
        if size(r1,2) ~= 1
            road_1 = [];
            for j = 1:size(r1,2)
                tube1 = r1(j);
                if ~any(ismember(predecessors_in{tube1},r1))
                    break
                end
            end
            road_1(1,end+1) = tube1;
            while size(road_1,2)~=size(r1,2)
                 a = find(ismember(predecessors_out{tube1},r1) == 1);
                road_1(1,end+1) = predecessors_out{tube1}(a);
                tube1 = predecessors_out{tube1}(a);
            end
        else
            road_1 = r1;
        end


        if size(r2,2) ~= 1
            road_2 = [];
            for j = 1:size(r2,2)
                tube2 = r2(j);
                if ~any(ismember(predecessors_in{tube2},r2))
                    road_2(1,end+1) = tube2;
                    break
                end
            end

            while size(road_2,2)~=size(r2,2)
                a = find(ismember(predecessors_out{tube2},r2) == 1);
                road_2(1,end+1) = predecessors_out{tube2}(a);
                tube2 = predecessors_out{tube2}(a);
            end
        else
            road_2 = r2;
        end

    pdropPaths = [pdropPaths,road_1,road_2];
    end
       

% ---------- 辅助函数：在给定管道集合中追踪串联路径 ----------
end