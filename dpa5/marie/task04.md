### Task 04:

```
function find_diameter(G = (V, E)):
  diameter = 0
  paralllel for each node u in V:
    max_dst = parallel_BFS(G, u)
    if max_dst > diameter:
      atomic diameter = max_dst
  return diameter
  
function parallel_BFS(G = (V, E), start):
    max_dst = 0
    
    distances = array of length n initialized to all -1
    current = empty node array
    next = empty node array
    
    push start to current
    distances[start] = 0
    
    while current is not empty:
      parallel for each u in current:
        for each v in N(u):
          if distances[v] == -1:
            atomic push v to next
            distances[v] = distances[u] + 1
            if distances[u] + 1 > max_dst:
              atomic max_dst = distances[u] + 1
      current = next
      next = empty node array
      
    return max_dst          
```

A parallel BFS algorithm runs in $O(nlogn)$ while using $O(n)$ processors. Since we can assume concurrent reads
we can run all n BFS instances at the same time using $O(n^2)$ processors for a total parallel runtime of $O(nlogn)$
as all BFSs all executed simultaneously.