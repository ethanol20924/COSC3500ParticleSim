# Particle Collision Simulation
*How to go zoom :) ~ Ethan Lo | s4638706*

## Research

So from milestone 1 we have a serial implementation of a particle simulation. It's really slow - obviously - so how do we make it go fast?

Milestone 2 requires us to implement parallel implementations of the simulation and profile it (I think, there's literally no task sheet what the fuck). So far the tools explored include:
- AVX
- OpenMP
- CUDA
- MPI

However, before we pick which one of these tools we wish to use, we must first figure out how we're going structure this simulation. In the serial implementation, particles were stored in just a flat vector. While this was nice as it puts everything in one contiguous block of memory, this will not suffice with a parallelised implementation. Due to the highly interconnected nature of the simulation (particles interacting with other particles), it is simply not enough to iterate over this flat structure. So what are our options?

## Quadtrees

The first structure that most people probably think of is a quadtree. By dynamically subdividing the simulation space into sub-spaces we can much more easily parallelise this by computing each leaf in parallel. The great thing about quadtrees as well is that it is very easy to find existing implementations:

[A simple and modern C++ quadtree implementation](https://github.com/pvigier/Quadtree)

This saves us the trouble of implementing a data structure from scratch (and since I'm not writing it, it's about 1000% more like to work). The author of this implmentation also has a [blog detailing his physics engine development](https://pvigier.github.io/2019/08/04/quadtree-collision-detection.html). Here he justifies his choice due to the versatility of the data structure, noting how it is not always the most optimal choice. [In another blog he references](https://0fps.net/2015/01/23/collision-detection-part-3-benchmarks/), various collision detection implementations are benchmarked. While trees are the fastest in most situations, grids are the fastest for uniformly distributed entities. This uniform distribution just happen to be exactly what the conditions of our simulation is.

## Grids

Ok so what about grids? Grids seem to be very popular, almost every post asking about quadtrees will get a response recommending they consider using grids (see [Efficient (and well explained) implementation of a Quadtree for 2D collision detection](https://stackoverflow.com/questions/41946007/efficient-and-well-explained-implementation-of-a-quadtree-for-2d-collision-det)). So what are they and how do they work? In an extensive answer to this post [What is a coarse and fine grid search?](https://stackoverflow.com/questions/59795569/what-is-a-coarse-and-fine-grid-search), describes the benefits of grids as supposed to tree based structures. He covers the use a multi-level grid system composed of a coarse and fine grid, where "the coarser grid contains the occupied cells of the finer grid, leading to two constant time queries". The finer grid is also a loose, allowing cell sizes to shrink and expand based on the size of the elements inserted.

The same answerer also describes a structure which uses a "medium-sized container per grid row". He motivates this by saying that single list structures could "lead to more heap allocations/deallocations", "translate to more cache misses" and "reduce the ability to process the grid in parallel". The following code exerpt describes such a structure:

```
struct GridRow
{
     struct Node
     {
         // Element data
         ...

         // Stores the index into the 'nodes' array of the next element
         // in the cell or the next removed element. -1 indicates end of
         // list.
         int next = -1;
     };

     // Stores all the nodes in the row.
     std::vector<Node> nodes;

     // Stores the index of the first removed element or -1
     // if there are no removed elements.
     int free_element = -1;
};
```

Such a structure allows for "constant-time insertions and removals while still allowing some degree of spatial locality".

Finally, he ends off the answer by introducing what he calls a "hierarchical bitset", as described in the following image:

![](images/Sparse%20Bitset.png)

## More on Grids

[Simple Fast Adaptive Grid to Accelerate Collision Detection between AABB of Particles](https://www.codeproject.com/Articles/5327631/Simple-Fast-Adaptive-Grid-to-Accelerate-Collision) also describes a grid based optimisation for particle collision. In this he describes two types of grid implementations: collision masks and center points.Collision masks store the element in all the cells where there could be collisions (based on the AABB of the element), while grids of center points only store the the element in the cell containing the center point. Storing the center points yields faster insertion and removal, while storing collision masks improves collision scanning. However, as long as the particles are small enough for a given grid size, the collision mask method becomes more effective.

![](images/Collision%20Mask%20Grid.png)

The article also outlines an adaptive grid implementation which subdivides cells which are too dense. However this addition only lead to a marginal increase in performance, with added complexity.

## So what do?

Basically my plan for the moment is to use a collision mask grid which is grouped by rows. We will also add a hierarchical bitset like structure for rows which don't contain intersecting AABBs. Who knows if this will speed it up at all, maybe for sparse sims.