#include <filesystem>
#include <cereal/archives/binary.hpp>
#include <fstream>
#include <globaltypes.hpp> 

rcs_store_t loadrcsstore(std::filesystem::path const & rcsstore_path)
{
    std::fstream rcsstream{rcsstore_path};
    if (!rcsstream.good())
      throw std::runtime_error{"Couldn`t read file."};

    rcs_store_t rcsstore{};
    {
        cereal::BinaryInputArchive rcsarchive{rcsstream};
        rcsstore.load(rcsarchive);
    }
    return rcsstore;
}