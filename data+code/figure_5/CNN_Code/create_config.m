function create_config(config_name, storage_location, image_list, number_pixels, offsets, classes, class_of_interest, original_images, file_prefix)
    config = struct();
    config.classes = classes;
    config.storage = storage_location;
    config.image_list = image_list;
    config.gauss_filt = 2;
    config.number_pixels = number_pixels;
    config.superpixels = struct();
    config.offsets = offsets;
    config.original_images = original_images;
    config.class_of_interest = class_of_interest;
    config.file_prefix = file_prefix;
    save([config_name, '_config'], 'config');
    disp('Structure Saved');
end