raw"""
    deform_mesh(mesh, displacement)
    deform_mesh(mesh, basis, amplitude)

Return a new [`Mesh`](@ref) with its vertices displaced, recomputing the panel
centres, normals, areas and radii from the moved geometry.

This is how shape sensitivities are taken on a hull that came out of a file.
The importers — `structured_surface_waterline_mesh`, the waterline clipping and
the STL reader — build their coordinates in `Float64`, so a dual number cannot
be pushed in through them. Deforming afterwards sidesteps that entirely: `Mesh`
is parameterised on its matrix types, so displacing an imported hull by
`amplitude * basis` with a dual `amplitude` yields a genuinely dual-valued mesh,
and `d/d(amplitude)` at zero is the shape derivative along `basis`.

The second form takes a fixed deformation field and a scalar amplitude, which
is the shape convenient for differentiation:

```julia
gradient = ForwardDiff.derivative(0.0) do amplitude
    volume(deform_mesh(hull, basis, amplitude))
end
```

`displacement` (or `basis`) has one row per *vertex*, matching `mesh.vertices`.
Panels whose corners become degenerate under the deformation keep their original
geometry rather than producing a zero-area panel, so a large amplitude
degrades rather than throwing; check `mesh_signed_volume` if that matters.

Vertices are not welded, so a field defined per vertex must give the same
displacement to every copy of a shared corner or the surface will tear. Build it
from `mesh.vertices` — any smooth function of position does this by
construction.
"""
function deform_mesh(mesh::Mesh, displacement::AbstractMatrix)
    size(displacement) == size(mesh.vertices) || throw(DimensionMismatch(
        "displacement must have size $(size(mesh.vertices))",
    ))
    vertices = mesh.vertices .+ displacement

    # Built without writing into preallocated arrays, so that this is
    # differentiable in REVERSE mode as well as forward. Reverse mode is the
    # useful direction for shape work — one adjoint pass gives the derivative
    # with respect to every vertex at once, where forward mode costs a pass per
    # variable — and it is the mode the published solver is built around.
    # Zygote cannot see through `setindex!`, so nothing here may mutate.
    panels = map(1:mesh.nfaces) do panel
        geometry = _quad_geometry(_panel_corners(vertices, mesh.faces, panel))
        isnothing(geometry) &&
            return (center = collect(@view mesh.centers[panel, :]),
                normal = collect(@view mesh.normals[panel, :]),
                area = mesh.areas[panel], radius = mesh.radii[panel])
        return (center = geometry.center, normal = geometry.normal,
            area = geometry.area, radius = geometry.radius)
    end

    centers = permutedims(reduce(hcat, map(panel -> panel.center, panels)))
    normals = permutedims(reduce(hcat, map(panel -> panel.normal, panels)))
    areas = map(panel -> panel.area, panels)
    radii = map(panel -> panel.radius, panels)
    return Mesh(vertices, mesh.faces, centers, normals, areas, radii,
        mesh.nvertices, mesh.nfaces)
end

# The four corners of a panel as a 4x3 matrix. Written with four explicit rows
# rather than a comprehension because Zygote has no adjoint for splatting a
# generator into `vcat`.
function _panel_corners(vertices, faces, panel::Integer)
    # Face indices are stored zero-based, as the importers write them.
    return vcat(transpose(@view vertices[faces[panel, 1] + 1, :]),
        transpose(@view vertices[faces[panel, 2] + 1, :]),
        transpose(@view vertices[faces[panel, 3] + 1, :]),
        transpose(@view vertices[faces[panel, 4] + 1, :]))
end

deform_mesh(mesh::Mesh, basis::AbstractMatrix, amplitude::Real) =
    deform_mesh(mesh, amplitude .* basis)

function deform_mesh(grid::StructuredPanelGrid, args...)
    return StructuredPanelGrid(deform_mesh(grid.mesh, args...), grid.strips)
end

raw"""
    beam_deformation_basis(mesh)

Deformation field that scales a hull's breadth, as a per-vertex displacement.

The simplest shape parameter worth a sensitivity, and the one whose derivative
can be checked against a hull generated parametrically at two breadths. Unit
amplitude doubles the breadth, so `deform_mesh(mesh, beam_deformation_basis(mesh),
a)` gives a hull of breadth ``(1+a)B``.
"""
function beam_deformation_basis(mesh::Mesh)
    basis = zero(mesh.vertices)
    basis[:, 2] .= @view mesh.vertices[:, 2]
    return basis
end

beam_deformation_basis(grid::StructuredPanelGrid) = beam_deformation_basis(grid.mesh)
