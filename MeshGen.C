#include "MeshGen.h"

using namespace std;

// Bring in everything from the libMesh namespace
using namespace libMesh;


MeshGen::MeshGen()
{

}

MeshGen::~MeshGen()
{

}

int MeshGen::idij(int nx, int i, int j)
{
    return i + j*(nx+1);
}

int MeshGen::idijk(int nx, int ny, int i, int j, int k)
{
    return i + (2*nx+1)*(j + k*(2*ny+1));
}


void MeshGen::create_mesh(double xmin, double xmax, double ymin, double ymax, int nx, int ny, Mesh & mesh)
{
    mesh.clear();
    mesh.reserve_nodes( (nx+1)*(ny+1) );

    BoundaryInfo & boundary_info = mesh.get_boundary_info();

    unsigned int node_id = 0;

    for (unsigned int j=0; j<=ny; j++)
        for (unsigned int i=0; i<=nx; i++)
            mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(nx),
                                  static_cast<Real>(j)/static_cast<Real>(ny),
                                  0.), node_id++);

    unsigned int elem_id = 0;
    for (unsigned int j=0; j<ny; j++)
        for (unsigned int i=0; i<nx; i++)
        {
            double x_cur = ((xmax-xmin)/(nx))*i+((xmax-xmin)/(nx))*0.5+xmin;

            // Add first Tri3
            Elem * elem = new Tri3;
            elem->set_id(elem_id++);
            elem = mesh.add_elem (elem);

            elem->set_node(0) = mesh.node_ptr(idij(nx,i,j)    );
            elem->set_node(1) = mesh.node_ptr(idij(nx,i+1,j)  );
            elem->set_node(2) = mesh.node_ptr(idij(nx,i+1,j+1));

            if (j == 0)
                boundary_info.add_side(elem, 0, 1);

            //if (i == (nx-1))
            //boundary_info.add_side(elem, 1, 1);

            // Add second Tri3
            elem = new Tri3;
            elem->set_id(elem_id++);
            elem = mesh.add_elem (elem);

            elem->set_node(0) = mesh.node_ptr(idij(nx,i,j)    );
            elem->set_node(1) = mesh.node_ptr(idij(nx,i+1,j+1));
            elem->set_node(2) = mesh.node_ptr(idij(nx,i,j+1)  );

            if (j == (ny-1))
            {
                if(fabs(x_cur)>=5.0)
                    boundary_info.add_side(elem, 1, 2);
                else
                    boundary_info.add_side(elem, 1, 3);
            }

            //if (i == 0)
            //boundary_info.add_side(elem, 2, 3);
        }

    for (unsigned int p=0; p<mesh.n_nodes(); p++)
    {
        mesh.node_ref(p)(0) = (mesh.node_ref(p)(0))*(xmax-xmin) + xmin;
        mesh.node_ref(p)(1) = (mesh.node_ref(p)(1))*(ymax-ymin) + ymin;
    }

    mesh.prepare_for_use();
}

void MeshGen::create_mesh_Cooks2D(double xmin, double xmax, double ymin, double ymax, int nx, int ny, Mesh & mesh)
{
    mesh.clear();
    mesh.reserve_nodes( (nx+1)*(ny+1) );

    BoundaryInfo & boundary_info = mesh.get_boundary_info();

    unsigned int node_id = 0;

    for (unsigned int j=0; j<=ny; j++)
        for (unsigned int i=0; i<=nx; i++)
            mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(nx),
                                  static_cast<Real>(j)/static_cast<Real>(ny),
                                  0.), node_id++);

    unsigned int elem_id = 0;
    for (unsigned int j=0; j<ny; j++)
        for (unsigned int i=0; i<nx; i++)
        {
            double x_cur = ((xmax-xmin)/(nx))*i+((xmax-xmin)/(nx))*0.5+xmin;

            // Add first Tri3
            Elem * elem = new Tri3;
            elem->set_id(elem_id++);
            elem = mesh.add_elem (elem);

            elem->set_node(0) = mesh.node_ptr(idij(nx,i,j)    );
            elem->set_node(1) = mesh.node_ptr(idij(nx,i+1,j)  );
            elem->set_node(2) = mesh.node_ptr(idij(nx,i+1,j+1));

            if (j == 0)
                boundary_info.add_side(elem, 0, 3);

            if (i == (nx-1))
                boundary_info.add_side(elem, 1, 2);

            // Add second Tri3
            elem = new Tri3;
            elem->set_id(elem_id++);
            elem = mesh.add_elem (elem);

            elem->set_node(0) = mesh.node_ptr(idij(nx,i,j)    );
            elem->set_node(1) = mesh.node_ptr(idij(nx,i+1,j+1));
            elem->set_node(2) = mesh.node_ptr(idij(nx,i,j+1)  );

            if (j == (ny-1))
            {
                //if(fabs(x_cur)>=5.0)
                boundary_info.add_side(elem, 1, 4);
                //else
                // boundary_info.add_side(elem, 1, 3);
            }

            if (i == 0)
                boundary_info.add_side(elem, 2, 1);
        }

    for (unsigned int p=0; p<mesh.n_nodes(); p++)
    {
        mesh.node_ref(p)(0) = (mesh.node_ref(p)(0))*(xmax-xmin) + xmin;
        mesh.node_ref(p)(1) = (mesh.node_ref(p)(1))*(ymax-ymin) + ymin;
    }

    for (unsigned int p=0; p<mesh.n_nodes(); p++)
    {

        mesh.node_ref(p)(1) = mesh.node_ref(p)(1)+(0.9167-0.13257*mesh.node_ref(p)(1))*mesh.node_ref(p)(0);
    }

    mesh.prepare_for_use();
}

void MeshGen::create_mesh_cube2D(double xmin, double xmax, double ymin, double ymax, int nx, int ny, Mesh & mesh)
{
    mesh.clear();
    mesh.reserve_nodes( (nx+1)*(ny+1) );

    BoundaryInfo & boundary_info = mesh.get_boundary_info();

    unsigned int node_id = 0;

    for (unsigned int j=0; j<=ny; j++)
        for (unsigned int i=0; i<=nx; i++)
            mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(nx),
                                  static_cast<Real>(j)/static_cast<Real>(ny),
                                  0.), node_id++);

    unsigned int elem_id = 0;
    for (unsigned int j=0; j<ny; j++)
        for (unsigned int i=0; i<nx; i++)
        {
            double x_cur = ((xmax-xmin)/(nx))*i+((xmax-xmin)/(nx))*0.5+xmin;

            // Add first Tri3
            Elem * elem = new Tri3;
            elem->set_id(elem_id++);
            elem = mesh.add_elem (elem);

            elem->set_node(0) = mesh.node_ptr(idij(nx,i,j)    );
            elem->set_node(1) = mesh.node_ptr(idij(nx,i+1,j)  );
            elem->set_node(2) = mesh.node_ptr(idij(nx,i+1,j+1));

            if (j == 0)
                boundary_info.add_side(elem, 0, 3);

            if (i == (nx-1))
                boundary_info.add_side(elem, 1, 2);

            // Add second Tri3
            elem = new Tri3;
            elem->set_id(elem_id++);
            elem = mesh.add_elem (elem);

            elem->set_node(0) = mesh.node_ptr(idij(nx,i,j)    );
            elem->set_node(1) = mesh.node_ptr(idij(nx,i+1,j+1));
            elem->set_node(2) = mesh.node_ptr(idij(nx,i,j+1)  );

            if (j == (ny-1))
            {
                //if(fabs(x_cur)>=5.0)
                boundary_info.add_side(elem, 1, 4);
                //else
                // boundary_info.add_side(elem, 1, 3);
            }

            if (i == 0)
                boundary_info.add_side(elem, 2, 1);
        }

    for (unsigned int p=0; p<mesh.n_nodes(); p++)
    {
        mesh.node_ref(p)(0) = (mesh.node_ref(p)(0))*(xmax-xmin) + xmin;
        mesh.node_ref(p)(1) = (mesh.node_ref(p)(1))*(ymax-ymin) + ymin;
    }

    mesh.prepare_for_use();
}

void MeshGen::create_mesh_cube_swing_2D(double xmin, double xmax, double ymin, double ymax, int nx, int ny, Mesh & mesh)
{
    mesh.clear();
    mesh.reserve_nodes( (nx+1)*(ny+1) );

    BoundaryInfo & boundary_info = mesh.get_boundary_info();

    unsigned int node_id = 0;

    for (unsigned int j=0; j<=ny; j++)
        for (unsigned int i=0; i<=nx; i++)
            mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(nx),
                                  static_cast<Real>(j)/static_cast<Real>(ny),
                                  0.), node_id++);

    unsigned int elem_id = 0;
    for (unsigned int j=0; j<ny; j++)
        for (unsigned int i=0; i<nx; i++)
        {
            double x_cur = ((xmax-xmin)/(nx))*i+((xmax-xmin)/(nx))*0.5+xmin;

            // Add first Tri3
            Elem * elem = new Tri3;
            elem->set_id(elem_id++);
            elem = mesh.add_elem (elem);

            elem->set_node(0) = mesh.node_ptr(idij(nx,i,j)    );
            elem->set_node(1) = mesh.node_ptr(idij(nx,i+1,j)  );
            elem->set_node(2) = mesh.node_ptr(idij(nx,i+1,j+1));

            if (j == 0)
                boundary_info.add_side(elem, 0, 1);

            if (i == (nx-1))
                boundary_info.add_side(elem, 1, 1);

            // Add second Tri3
            elem = new Tri3;
            elem->set_id(elem_id++);
            elem = mesh.add_elem (elem);

            elem->set_node(0) = mesh.node_ptr(idij(nx,i,j)    );
            elem->set_node(1) = mesh.node_ptr(idij(nx,i+1,j+1));
            elem->set_node(2) = mesh.node_ptr(idij(nx,i,j+1)  );

            if (j == (ny-1))
            {
                //if(fabs(x_cur)>=5.0)
                boundary_info.add_side(elem, 1, 1);
                //else
                // boundary_info.add_side(elem, 1, 3);
            }

            if (i == 0)
                boundary_info.add_side(elem, 2, 1);
        }

    for (unsigned int p=0; p<mesh.n_nodes(); p++)
    {
        mesh.node_ref(p)(0) = (mesh.node_ref(p)(0))*(xmax-xmin) + xmin;
        mesh.node_ref(p)(1) = (mesh.node_ref(p)(1))*(ymax-ymin) + ymin;
    }

    mesh.prepare_for_use();
}

void MeshGen::create_mesh_Cooks3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int nx, int ny, int nz, Mesh & mesh)
{
    mesh.clear();
    mesh.reserve_elem(nx*ny*nz);
    mesh.reserve_nodes( (2*nx+1)*(2*ny+1)*(2*nz+1) );

    BoundaryInfo & boundary_info = mesh.get_boundary_info();

    unsigned int node_id = 0;

    for (unsigned int k=0; k<=(2*nz); k++)
        for (unsigned int j=0; j<=(2*ny); j++)
            for (unsigned int i=0; i<=(2*nx); i++)
                mesh.add_point(Point(static_cast<Real>(i)/static_cast<Real>(2*nx),
                                     static_cast<Real>(j)/static_cast<Real>(2*ny),
                                     static_cast<Real>(k)/static_cast<Real>(2*nz)), node_id++);

    unsigned int elem_id = 0;

    for (unsigned int k=0; k<(2*nz); k += 2)
        for (unsigned int j=0; j<(2*ny); j += 2)
            for (unsigned int i=0; i<(2*nx); i += 2)
            {
                Elem * elem = static_cast<Elem *>(new Hex27);
                elem->set_id(elem_id++);
                elem = mesh.add_elem (elem);

                elem->set_node(0)  = mesh.node_ptr(idijk(nx,ny,i,  j,  k)  );
                elem->set_node(1)  = mesh.node_ptr(idijk(nx,ny,i+2,j,  k)  );
                elem->set_node(2)  = mesh.node_ptr(idijk(nx,ny,i+2,j+2,k)  );
                elem->set_node(3)  = mesh.node_ptr(idijk(nx,ny,i,  j+2,k)  );
                elem->set_node(4)  = mesh.node_ptr(idijk(nx,ny,i,  j,  k+2));
                elem->set_node(5)  = mesh.node_ptr(idijk(nx,ny,i+2,j,  k+2));
                elem->set_node(6)  = mesh.node_ptr(idijk(nx,ny,i+2,j+2,k+2));
                elem->set_node(7)  = mesh.node_ptr(idijk(nx,ny,i,  j+2,k+2));
                elem->set_node(8)  = mesh.node_ptr(idijk(nx,ny,i+1,j,  k)  );
                elem->set_node(9)  = mesh.node_ptr(idijk(nx,ny,i+2,j+1,k)  );
                elem->set_node(10) = mesh.node_ptr(idijk(nx,ny,i+1,j+2,k)  );
                elem->set_node(11) = mesh.node_ptr(idijk(nx,ny,i,  j+1,k)  );
                elem->set_node(12) = mesh.node_ptr(idijk(nx,ny,i,  j,  k+1));
                elem->set_node(13) = mesh.node_ptr(idijk(nx,ny,i+2,j,  k+1));
                elem->set_node(14) = mesh.node_ptr(idijk(nx,ny,i+2,j+2,k+1));
                elem->set_node(15) = mesh.node_ptr(idijk(nx,ny,i,  j+2,k+1));
                elem->set_node(16) = mesh.node_ptr(idijk(nx,ny,i+1,j,  k+2));
                elem->set_node(17) = mesh.node_ptr(idijk(nx,ny,i+2,j+1,k+2));
                elem->set_node(18) = mesh.node_ptr(idijk(nx,ny,i+1,j+2,k+2));
                elem->set_node(19) = mesh.node_ptr(idijk(nx,ny,i,  j+1,k+2));

                elem->set_node(20) = mesh.node_ptr(idijk(nx,ny,i+1,j+1,k)  );
                elem->set_node(21) = mesh.node_ptr(idijk(nx,ny,i+1,j,  k+1));
                elem->set_node(22) = mesh.node_ptr(idijk(nx,ny,i+2,j+1,k+1));
                elem->set_node(23) = mesh.node_ptr(idijk(nx,ny,i+1,j+2,k+1));
                elem->set_node(24) = mesh.node_ptr(idijk(nx,ny,i,  j+1,k+1));
                elem->set_node(25) = mesh.node_ptr(idijk(nx,ny,i+1,j+1,k+2));
                elem->set_node(26) = mesh.node_ptr(idijk(nx,ny,i+1,j+1,k+1));


                if (k == 0)
                    boundary_info.add_side(elem, 0, 0);

                if (k == 2*(nz-1))
                    boundary_info.add_side(elem, 5, 5);

                if (j == 0)
                    boundary_info.add_side(elem, 1, 1);

                if (j == 2*(ny-1))
                    boundary_info.add_side(elem, 3, 3);

                if (i == 0)
                    boundary_info.add_side(elem, 4, 4);

                if (i == 2*(nx-1))
                    boundary_info.add_side(elem, 2, 2);
            }


    std::vector<Elem *> new_elements;
    new_elements.reserve(24*mesh.n_elem());

    MeshBase::element_iterator       el     = mesh.elements_begin();
    const MeshBase::element_iterator end_el = mesh.elements_end();

    for ( ; el != end_el;  ++el)
    {
        // Get a pointer to the HEX27 element.
        Elem * base_hex = *el;

        // Get a pointer to the node located at the HEX27 centroid
        Node * apex_node = base_hex->node_ptr(26);

        // Container to catch ids handed back from BoundaryInfo
        std::vector<boundary_id_type> ids;

        for (unsigned short s=0; s<base_hex->n_sides(); ++s)
        {
            // Get the boundary ID(s) for this side
            boundary_info.boundary_ids(*el, s, ids);

            // We're creating this Mesh, so there should be 0 or 1 boundary IDs.
            libmesh_assert(ids.size() <= 1);

            // A convenient name for the side's ID.
            boundary_id_type b_id = ids.empty() ? BoundaryInfo::invalid_id : ids[0];

            // Need to build the full-ordered side!
            UniquePtr<Elem> side = base_hex->build_side_ptr(s);


            // Build 4 sub-tets per side
            for (unsigned int sub_tet=0; sub_tet<4; ++sub_tet)
            {
                new_elements.push_back( new Tet4 );
                Elem * sub_elem = new_elements.back();
                sub_elem->set_node(0) = side->node_ptr(sub_tet);
                sub_elem->set_node(1) = side->node_ptr(8);                           // centroid of the face
                sub_elem->set_node(2) = side->node_ptr(sub_tet==3 ? 0 : sub_tet+1 ); // wrap-around
                sub_elem->set_node(3) = apex_node;                                   // apex node always used!

                // If the original hex was a boundary hex, add the new sub_tet's side
                // 0 with the same b_id.  Note: the tets are all aligned so that their
                // side 0 is on the boundary.
                if (b_id != BoundaryInfo::invalid_id)
                    boundary_info.add_side(sub_elem, 0, b_id);
            }

        }
    }

	{
		MeshBase::element_iterator       el     = mesh.elements_begin();
		const MeshBase::element_iterator end_el = mesh.elements_end();

		for ( ; el != end_el;  ++el)
		{
			boundary_info.remove(*el); // Safe even if *el has no boundary info.
			mesh.delete_elem(*el);
		}
	}

    for (std::size_t i=0; i<new_elements.size(); ++i)
        mesh.add_elem(new_elements[i]);
        
    for (unsigned int p=0; p<mesh.n_nodes(); p++)
    {
        mesh.node_ref(p)(0) = (mesh.node_ref(p)(0))*(xmax-xmin) + xmin;
        mesh.node_ref(p)(1) = (mesh.node_ref(p)(1))*(ymax-ymin) + ymin;
    }

    for (unsigned int p=0; p<mesh.n_nodes(); p++)
    {

        mesh.node_ref(p)(1) = mesh.node_ref(p)(1)+(0.9167-0.13257*mesh.node_ref(p)(1))*mesh.node_ref(p)(0);
    }

    mesh.prepare_for_use();
}

