#pragma once
#include "allscale/api/user/data/grid.h"
#include "amdados/app/static_grid.h"

//using namespace allscale::api::user;  // Albert commented

namespace amdados {
namespace app {


    template<unsigned ... sizes>
    struct layer;

     template<typename ... Layers>
     struct layers {
         enum {
             num_layers = sizeof...(Layers)
         };
     };

    template<unsigned ... sizes>
    struct size;

    // structures for each Cell configuration and number of layers for nesting (limited to 3)
    template<typename Layers>
    struct CellConfig {
        enum {
            num_layers = Layers::num_layers
        };
    };

    enum Direction {
        Up, Down, Left, Right
    };


    template<typename T, typename Size, typename Layers>
    struct GridLayerData;

    template<typename T, unsigned ... Sizes, unsigned ... Dims, typename ... Rest>
    struct GridLayerData<T, size<Sizes...>, layers<layer<Dims...>,Rest...> > {
        using data_type = allscale::utils::grid<T,Sizes...>;
        using nested_type = GridLayerData<T,size<Sizes*Dims...>,layers<Rest...>>;

        enum { layer_number = sizeof...(Rest)+1 };

        // the values to be stored on this layer
        data_type data;

        // the nested layers
        nested_type nested;

        unsigned getLayerNumber() const {
            return layer_number;
        }

        template<unsigned Layer>
        typename std::enable_if<Layer == layer_number, data_type&>::type
        getLayer() {
            return data;
        }

        template<unsigned Layer>
        typename std::enable_if<Layer == layer_number, const data_type&>::type
        getLayer() const {
            return data;
        }

        template<unsigned Layer>
        typename std::enable_if<Layer < layer_number,
            decltype(nested.template getLayer<Layer>())>::type
        getLayer() {
            return nested.template getLayer<Layer>();
        }

        template<unsigned Layer>
        typename std::enable_if<Layer < layer_number,
            decltype(static_cast<const nested_type&>(nested).template getLayer<Layer>())>::type
        getLayer() const {
            return nested.template getLayer<Layer>();
        }

        template<typename Op>
        void forAllOnLayer(unsigned layer, const Op& op) {
            if (layer == getLayerNumber()) {
                // apply it to this value)
                data.for_each(op);
            } else {
                nested.forAllOnLayer(layer,op);
            }
        }

        // returns vector of boundary data in each direction
        // TODO: avoid the temporary array creation.
        std::vector<T> getBoundary(unsigned layer, Direction dir) const {
            if (layer == getLayerNumber()) {
                int size[] = {Sizes...};
                int Nx = size[0];
                int Ny = size[1];
                switch (dir) {
                    case Up: {      // returns data from TOP strip of domain to neighbour
                        std::vector<T> res(Nx);
                        for (int x = 0; x < Nx; ++x) res[x] = data[{x,Ny-1}];
                        return res;
                    }
                    case Down: {    // returns data from BOTTOM strip of domain to neighbour
                        std::vector<T> res(Nx);
                        for (int x = 0; x < Nx; ++x) res[x] = data[{x,0}];
                        return res;
                    }
                    case Left: {
                        std::vector<T> res(Ny);
                        for (int y = 0; y < Ny; ++y) res[y] = data[{0,y}];
                        return res;
                    }
                    case Right: {
                        std::vector<T> res(Ny);
                        for (int y = 0; y < Ny; ++y) res[y] = data[{Nx-1,y}];
                        return res;
                    }
                }
            }
            return nested.getBoundary(layer, dir);
        }

        void setBoundary(unsigned layer, Direction dir, const std::vector<T> & boundary) {
            if (layer == getLayerNumber()) {
                int size[] = {Sizes...};
                int Nx = size[0];
                int Ny = size[1];
                switch (dir) {
                    case Up: {
                        assert_true(boundary.size() == static_cast<size_t>(Nx));
                        for (int x = 0; x < Nx; ++x) data[{x,Ny-1}] = boundary[x];
                        return;
                    }
                    case Down: {
                        assert_true(boundary.size() == static_cast<size_t>(Nx));
                        for (int x = 0; x < Nx; ++x) data[{x,0}] = boundary[x];
                        return;
                    }
                    case Left: {
                        assert_true(boundary.size() == static_cast<size_t>(Ny));
                        for (int y = 0; y < Ny; ++y) data[{0,y}] = boundary[y];
                        return;
                    }
                    case Right: {
                        assert_true(boundary.size() == static_cast<size_t>(Ny));
                        for (int y = 0; y < Ny; ++y) data[{Nx-1,y}] = boundary[y];
                        return;
                    }
                }
            }
            nested.setBoundary(layer, dir, boundary);
        }
    };





    template<typename T, typename CellConfig>
    struct Cell;


    template<typename T, typename Layers>
    struct Cell<T,CellConfig<Layers>> {

        // determines the active layer of this grid cell
        unsigned active_layer = 0;

        // the data stored in
        GridLayerData<T,size<1,1>,Layers> data;

        Cell& operator=(const Cell& other) {
            if(this == &other) return *this;
            active_layer = other.active_layer;
            switch(active_layer) {
            case 0: getLayer<0>() = other.getLayer<0>(); return *this;
            case 1: getLayer<1>() = other.getLayer<1>(); return *this;
            case 2: getLayer<2>() = other.getLayer<2>(); return *this;
            case 3: getLayer<3>() = other.getLayer<3>(); return *this;
            default: assert_true(false && "Unsupported number of levels!");
            }
            return *this;
        }

        void setActiveLayer(unsigned level) {
            active_layer = level;
        }



        unsigned getActiveLayer() const {
            return active_layer;
        }

        template<unsigned Layer>
        auto getLayer() -> decltype(data.template getLayer<Layer>())& {
            return data.template getLayer<Layer>();
        }

        template<unsigned Layer>
        auto getLayer() const -> const decltype(data.template getLayer<Layer>())& {
            return data.template getLayer<Layer>();
        }

        template<typename Op>
        void forAllActiveNodes(const Op& op) {
            data.forAllOnLayer(active_layer, op);
            //data.DiscretizeElements(active_layer);
        }


        std::vector<T> getBoundary(Direction dir) const {
            return data.getBoundary(active_layer,dir);
        }

        void setBoundary(Direction dir, const std::vector<T>& boundary) {
            data.setBoundary(active_layer,dir,boundary);
        }



//      template<unsigned Layer>
//      std::enable_if<Layers::num_layers -
//      T& get(const utils::Coordinate<3>& pos) {
//
//      }

    };


    template<typename T, typename Size, typename Layers>
    struct GridLayerData;

    template<typename T,unsigned ... Sizes>
    struct GridLayerData<T,size<Sizes...>,layers<>> {

        using data_type =  allscale::utils::grid<T,Sizes...>;

        // the values to be stored on this last layer
        data_type data;

        unsigned getLayerNumber() const {
            return 0;
        }

        template<unsigned Layer>
        typename std::enable_if<Layer == 0, data_type&>::type
        getLayer() {
            return data;
        }

        template<unsigned Layer>
        typename std::enable_if<Layer == 0, const data_type&>::type
        getLayer() const {
            return data;
        }

        template<typename Op>
        void forAllOnLayer(unsigned layer, const Op& op) {
            assert_eq(layer, 0) << "Error: trying to access layer " << layer << " --no such layer!";
            data.for_each(op);
        }

        // TODO: avoid the temporary array creation.
        std::vector<T> getBoundary(unsigned layer, Direction dir) const {
            if (layer == 0) {
                int size[] = {Sizes...};
                int Nx = size[0];
                int Ny = size[1];
                switch (dir) {
                    case Up: {      // returns data from TOP strip of domain to neighbour
                        std::vector<T> res(Nx);
                        for (int x = 0; x < Nx; ++x) res[x] = data[{x,Ny-1}];
                        return res;
                    }
                    case Down: {    // returns data from BOTTOM strip of domain to neighbour
                        std::vector<T> res(Nx);
                        for (int x = 0; x < Nx; ++x) res[x] = data[{x,0}];
                        return res;
                    }
                    case Left: {
                        std::vector<T> res(Ny);
                        for (int y = 0; y < Ny; ++y) res[y] = data[{0,y}];
                        return res;
                    }
                    case Right: {
                        std::vector<T> res(Ny);
                        for (int y = 0; y < Ny; ++y) res[y] = data[{Nx-1,y}];
                        return res;
                    }
                }
            }
            assert_true(false && "No such layer!");
            return std::vector<T>();
        }

        void setBoundary(unsigned layer, Direction dir, const std::vector<T> & boundary) {
            if (layer == 0) {
                int size[] = {Sizes...};
                int Nx = size[0];
                int Ny = size[1];
                switch (dir) {
                    case Up: {
                        assert_true(boundary.size() == static_cast<size_t>(Nx));
                        for (int x = 0; x < Nx; ++x) data[{x,Ny-1}] = boundary[x];
                        return;
                    }
                    case Down: {
                        assert_true(boundary.size() == static_cast<size_t>(Nx));
                        for (int x = 0; x < Nx; ++x) data[{x,0}] = boundary[x];
                        return;
                    }
                    case Left: {
                        assert_true(boundary.size() == static_cast<size_t>(Ny));
                        for (int y = 0; y < Ny; ++y) data[{0,y}] = boundary[y];
                        return;
                    }
                    case Right: {
                        assert_true(boundary.size() == static_cast<size_t>(Ny));
                        for (int y = 0; y < Ny; ++y) data[{Nx-1,y}] = boundary[y];
                        return;
                    }
                }
            }
            assert_true(false && "No such layer!");
        }
    };

} // end namespace app
} // end namespace amdados

