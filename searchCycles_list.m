function cycleList = searchCycles(edgeMap)

    tic
    global graph cycles numCycles cycles_dist;
    graph = edgeMap;
    numCycles = 0;
    cycles = {};
    cycles_dist = zeros(100000,1);
    for i = 1:size(graph,1)
        for j = 1:2
            findNewCycles(graph(i,j))
        end
    end
    % print out all found cycles
    for i = 1:size(cycles_dist,1)
        cycles_dist(i)
    end
    % return the result
    cycleList = cycles_dist;
    xlswrite('result1.xlsx',cycles_dist)
    toc

function findNewCycles(path)

global graph cycles numCycles cycles_dist;
startNode = path(1);
nextNode = nan;
sub = [];

% visit each edge and each node of each edge
for i = 1:size(graph,1)
    node1 = graph(i,1);
    node2 = graph(i,2);
       if (node1 == startNode) || (node2==startNode) %% this if is required
        if node1 == startNode
            nextNode = node2;
        elseif node2 == startNode
            nextNode = node1;
        end
        if ~(visited(nextNode, path))
            % neighbor node not on path yet
            sub = nextNode;
            sub = [sub path];
            % explore extended path
            findNewCycles(sub);
        elseif size(path,2) > 2 && nextNode == path(end)
            % cycle found
            p = rotate_to_smallest(path);
            inv = invert(p);
            if isNew(p) && isNew(inv)
                numCycles = numCycles + 1;
                cycles{numCycles} = p;
                cycles_dist(size(path,2))=cycles_dist(size(path,2))+1;
            end
        end
    end
end

function inv = invert(path)
    inv = rotate_to_smallest(path(end:-1:1));

% rotate cycle path such that it begins with the smallest node
function new_path = rotate_to_smallest(path)
    [~,n] = min(path);
    new_path = [path(n:end), path(1:n-1)];

function result = isNew(path)
    global cycles
    result = 1;
    for i = 1:size(cycles,2)
        if size(path,2) == size(cycles{i},2) && all(path == cycles{i})
            result = 0;
            break;
        end
    end

function result = visited(node,path)
    result = 0;
    if isnan(node) && any(isnan(path))
        result = 1;
        return
    end
    for i = 1:size(path,2)
        if node == path(i)
            result = 1;
            break
        end
    end