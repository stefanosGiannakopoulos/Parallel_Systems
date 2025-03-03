#include <stdio.h>
#include <unistd.h>

int main() {
    // Get the page size using sysconf
    long page_size = sysconf(_SC_PAGESIZE);

    if (page_size == -1) {
        perror("sysconf");
        return 1;
    }

    // Print the page size
    printf("Page size: %ld bytes\n", page_size);
    
    return 0;
}

