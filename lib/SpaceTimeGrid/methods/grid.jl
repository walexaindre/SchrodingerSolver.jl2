export get_indexed_elements, evaluate_indexed_elements

function linear_indexing(Grid::SpaceTimeGrid)
    meta = get_metadata(Grid)
    max_size = linear_size(meta)

    linear_index_wrap(index) = linear_indexing(index, meta)

    indexes = collect(1:max_size)

    linear_index_wrap.(indexes)

end

function get_indexed_elements(Grid::SpaceTimeGrid)
    function get_element_wrap(index)
        get_element(Grid, index...)
    end
    indexing = linear_indexing(Grid)
    get_element_wrap.(indexing)
end

function evaluate_indexed_elements(Grid::SpaceTimeGrid, fun::Function)
    indexed_collection = get_indexed_elements(Grid)

    function evaluate_wrap(tuple)
        fun(tuple...)
    end

    evaluate_wrap.(indexed_collection)
end